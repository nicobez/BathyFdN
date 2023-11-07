########################################
# SPDE Non Stationary kriging (KO? KS?) with the ML parameters 
# myOptimNoStatParamInit (or myOptimNoStatParam if a full optimization has been done)
########################################

sel <- myGrid['delatedPolygon']
selRank <- (1:length(sel))[sel==1]
temp <- rep(NA,length(sel))

myModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=c(1,1),sill=1)

propNugget <- myOptimNoStatParamInit[1]
r1Start <- myOptimNoStatParamInit[2]
r1End <- myOptimNoStatParamInit[3]
r2Start <- myOptimNoStatParamInit[4]
r2End <- myOptimNoStatParamInit[5]

myGrid['r1'] <- setAnisoRange(r1Start,r1End) 
myGrid['r2'] <- setAnisoRange(r2Start,r2End) 
myGrid['sill'] <- driftSd(myGrid['distAcross'])**2
myGrid$setLocator("dirAcross",ELoc_NOSTAT(),0)
myGrid$setLocator("r1",ELoc_NOSTAT(),1)
myGrid$setLocator("r2",ELoc_NOSTAT(),2)
#myGrid$setLocator("sill",ELoc_NOSTAT(),3)

#nostat = NoStatArray(c("A1","R1","R2","V1"),myGrid)
nostat = NoStatArray(c("A1","R1","R2"),myGrid)
Model_addNoStat(myModel,nostat)

myDb['sill'] <- driftSd(myDb['distAcross'])**2
myDb['nugget'] <- propNugget * myDb['sill']

Q = cs_toTL(PrecisionOpCs(myMesh,myModel)$getQ()) #Precision matrix for sill=1
Dsill <- Diagonal(length(selRank),sqrt(myGrid['sill'][selRank]))
Q <- Dsill %*% Q %*% Dsill # precision matrix for non stationary sills

# if nugget is a scalar
# M1 <- nugget^2*Q + t(Aproj)%*%Aproj
# M2 <- t(Aproj)%*%myDb["depthC"] #myDb["depthResiduals"]

# if the nugget is not a scalar
D <- Diagonal(length(myDb['nugget']),myDb['nugget'])
invD <- Diagonal(length(myDb['nugget']),1/myDb['nugget'])
M1 <- Q + t(Aproj) %*% invD %*% Aproj
M2 <- t(Aproj) %*% invD %*% myDb["depthC"] 
zKri <-  Matrix:::solve(M1,M2)

# affecting the result to the grid
temp[selRank] <-zKri
myGrid['kriNonStatML'] <- temp + driftM(myGrid['distAcross'])

f.image <- function(grid=myGrid,varName,col=rev(myPalette),...){
  image(unique(grid['x1']), unique(grid['x2']),
        matrix(ifelse(grid['Polygon']==1,grid[varName],NA),nx),
        asp=1,col=col,
        xlab="Longitude",ylab="Latitude",las=1,
        xlim=range(grid['x1'][grid['Polygon']==1]),...)
  lines(ref.line)
}

f.image(varName='kriNonStatML',zlim=c(0,75))

# myGrid$setLocator("Polygon",ELoc_SEL(),0)
# p = ggplot() + plot.grid(myGrid, name_raster='kriNonStatML',palette='viridis')
# ggPrint(p)
# myGrid$setLocator("delatedPolygon",ELoc_SEL(),0)

dev.print(device = png, file = "Res/krigingSPDEwithML.png",width=450,height=450)

###########################################
# Kriging estimation variance Using NonCond simulations
###########################################
cholQ = chol(Q)

nSimu <- 400 # 40 secondes for nSimu = 200
simuCond <- array(NA,dim=c(myGrid$getActiveSampleNumber(),nSimu))
for(i in 1:nSimu){
  #simuNonCondVertices = solve(cholQ,rnorm(NROW(Q)))[,1]
  simuNonCondVertices = solve(cholQ,rnorm(NROW(Q),0,sqrt(myGrid['sill'][selRank])))[,1]
  simuNonCondData <- Aproj%*%simuNonCondVertices + rnorm(NROW(Aproj),0,sqrt(myDb['nugget']))
  #M2 <- t(Aproj)%*%simuNonCondData
  M2 <- t(Aproj)%*%invD%*%simuNonCondData 
  
  nonCondKri <- Matrix:::solve(M1,M2)
  simuCondVertices <- zKri + simuNonCondVertices - nonCondKri
  simuCond[,i] <- as.vector(simuCondVertices)
}

simuCond <- simuCond + driftM(myGrid['distAcross'][selRank])
simuCond[simuCond < min(myDb['depth'],na.rm=T)] <- min(myDb['depth'],na.rm=T)

estimKriDepth <- apply(simuCond,1,"mean")
sdKriDepth <- apply(simuCond,1,"sd")

temp[selRank] <- sdKriDepth
sdKriDepth <- temp # * driftSd(myGrid['distAcross'])

temp[selRank] <- estimKriDepth
estimKriDepth <- temp #+ driftM(myGrid['distAcross'])
estimKriDepth[estimKriDepth < min(myDb['depth'],na.rm=T)] <- min(myDb['depth'],na.rm=T)

for(i in 1:10){
  temp[selRank] <- simuCond[,i]
  #temp <- temp + driftM(myGrid['distAcross'])
  #temp[temp < min(myDb['depth'],na.rm=T)] <- min(myDb['depth'],na.rm=T)
  assign(paste0("simuCondDepth",i),temp)
}

#myGrid['residualKriDepthByCondSimu'] <- tempKri
myGrid['estimKriDepthByCondSimu'] <- estimKriDepth
myGrid['sdKriDepthByCondSimu'] <- sdKriDepth
myGrid['cvKriDepthByCondSimu'] <- sdKriDepth/estimKriDepth
for(i in 1:10) myGrid[paste0("simuCondDepth",i)] <- get(paste0("simuCondDepth",i))

f.image(varName='estimKriDepthByCondSimu',zlim=c(0,75))
f.image(varName='sdKriDepthByCondSimu',col=myPalette, zlim=c(0,5))
f.image(varName='cvKriDepthByCondSimu',zlim=c(0,0.2))
f.image(varName='simuCondDepth1',zlim=c(0,100))
f.image(varName='simuCondDepth2',zlim=c(0,100))
f.image(varName='simuCondDepth3',zlim=c(0,100))
f.image(varName='simuCondDepth4',zlim=c(0,100))

image(unique(myGrid['x1']), unique(myGrid['x2']),
      matrix(ifelse(myGrid['Polygon']==1 & abs(myGrid['cvKriDepthByCondSimu']) < 0.05,myGrid['estimKriDepthByCondSimu'],NA),nx),
      zlim=c(0,200),
      asp=1,col=myPalette,
      xlab="Longitude",ylab="Latitude",las=1,
      xlim=range(myGrid['x1'][myGrid['Polygon']==1]))
lines(ref.line)

tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("simuCondDepth",1)],NA),nx))
plot(-tempRes$z[140,],ylim=c(-100,0),type="l")
for(i in 2:10){
  tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("simuCondDepth",i)],NA),nx))
lines(-tempRes$z[140,],col=i)
}
                

#### Polygon selection

myGridExtPolygon <- myGridExt$clone()
db_polygon(myGridExtPolygon, myPolygon, flag_nested=TRUE)

p = ggplot() + plot.grid(myGridExtPolygon, name_raster='pseudoKriDepth',palette='viridis')
ggPrint(p)

tempPolygon <- myGridExtPolygon['Polygon']
tempKri <- myGridExtPolygon['pseudoKriDepth']
tempKri[tempPolygon==0] <- NA
tempSd <- myGridExtPolygon['sdKriDepth']
tempSd[tempPolygon==0] <- NA
  
image(unique(myGridExt['x1']), unique(myGridExt['x2']),matrix(tempKri,275),asp=1,col=myPalette,
      xlab="Longitude",ylab="Latitude",las=1,
      xlim=range(myGridExtPolygon['x1'][tempPolygon==1]))
lines(ref.line)
dev.print(device = png, file = "Res/krigingSPDE.png",width=450,height=450)

image(unique(myGridExt['x1']), unique(myGridExt['x2']),matrix(tempSd,275),asp=1,col=myPalette,
      xlab="Longitude",ylab="Latitude",las=1,
      xlim=range(myGridExtPolygon['x1'][tempPolygon==1]))
dev.print(device = png, file = "Res/sdSPDE.png",width=450,height=450)


#######

res <- array(NA,dim=c(1000,5))
for(i in 1:1000) res[i,] <- rnorm(5,0,1:5)
hist(res[,1])
hist(res[,5])
