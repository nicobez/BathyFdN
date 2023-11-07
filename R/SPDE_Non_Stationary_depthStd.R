########################################
# SPDE Non Stationary kriging (KO? KS?) with the ML parameters myOptimNoStatParamStd
########################################
ndat <- length(myDb['depthStd'])

sel <- myGrid['delatedPolygon']
selRank <- (1:length(sel))[sel==1]
temp <- rep(NA,length(sel))

myModelNoStat = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=c(1,1),sill=1)

nugStart <- myOptimNoStatParamStd[1]
nugEnd <- myOptimNoStatParamStd[2]
r1Start <- myOptimNoStatParamStd[3]
r1End <- myOptimNoStatParamStd[4]
r2Start <- myOptimNoStatParamStd[5]
r2End <- myOptimNoStatParamStd[6]

myGrid['r1'] <- setAnisoRange(r1Start,r1End) 
#myGrid['r2'] <- setAnisoRange(r2Start,r2End) 
myGrid['r2'] <- setAnisoRange2(myGrid['distAcross']) 
myGrid['sill'] <- driftSd(myGrid['distAcross'])**2
myGrid$setLocator("dirAcross",ELoc_NOSTAT(),0)
myGrid$setLocator("r1",ELoc_NOSTAT(),1)
myGrid$setLocator("r2",ELoc_NOSTAT(),2)
#myGrid$setLocator("sill",ELoc_NOSTAT(),3)

#nostat = NoStatArray(c("A1","R1","R2","V1"),myGrid)
nostat = NoStatArray(c("A1","R1","R2"),myGrid)
Model_addNoStat(myModelNoStat,nostat)

QNoStat = cs_toTL(PrecisionOpCs(myMesh,myModelNoStat)$getQ()) #Precision matrix for sill=1

distAcross <- myDb['distAcross']
myDb['nugget'] <- (myDb['distAcross'] < 0) * 0 +
  (myDb['distAcross'] >= 0 & myDb['distAcross'] < 1) * (nugStart + (nugEnd-nugStart)*myDb['distAcross']) +
  (myDb['distAcross'] > 1) * nugEnd
DNoStat <- Diagonal(length(myDb['nugget']),myDb['nugget'])
invDNoStat <- Diagonal(length(myDb['nugget']),1/myDb['nugget'])

M1 <- QNoStat + t(A) %*% invDNoStat %*% A
M2 <- t(A) %*% invDNoStat %*% myDb["depthStd"] 
zKri <-  Matrix:::solve(M1,M2)

# affecting the result to the grid
temp[selRank] <-zKri
myGrid['nonStatKri'] <- driftSd(myGrid['distAcross'])*temp + driftM(myGrid['distAcross'])

###########################################
# Kriging estimation variance Using NonCond simulations
###########################################
cholQNoStat = chol(QNoStat)

nSimu <- 400 # 40 secondes for nSimu = 200
simuCond <- array(NA,dim=c(myGrid$getActiveSampleNumber(),nSimu))
for(i in 1:nSimu){
  simuNonCondVertices = solve(cholQNoStat,rnorm(NROW(QNoStat),0,1))[,1]
  simuNonCondData <- A %*% simuNonCondVertices + rnorm(NROW(A),0,sqrt(myDb['nugget']))
  M2 <- t(A) %*% invDNoStat %*% simuNonCondData 
  nonCondKri <- Matrix:::solve(M1,M2)
  simuCondVertices <- zKri + simuNonCondVertices - nonCondKri
  simuCond[,i] <- as.vector(simuCondVertices)
}

# from Y to Z= Y*sd + m
simuCond <- driftSd(myGrid['distAcross'][selRank])*simuCond + driftM(myGrid['distAcross'][selRank])
# truncation of Z to the observed extremum
simuCond[simuCond < min(myDb['depth'],na.rm=T)] <- min(myDb['depth'],na.rm=T)

# Affecting outputs to the grid
sdKriDepth <- apply(simuCond,1,"sd")
temp[selRank] <- sdKriDepth
sdKriDepth <- temp

for(i in 1:10){
  temp[selRank] <- simuCond[,i]
  assign(paste0("nonStatSimu",i),temp)
}

myGrid['nonStatKriSD'] <- sdKriDepth
myGrid['nonStatKriCV'] <- sdKriDepth/myGrid['nonStatKri']
for(i in 1:10) myGrid[paste0("nonStatSimu",i)] <- get(paste0("nonStatSimu",i))

f.image(varName='nonStatKri',polyName='myPoly',zlim=c(0,200)) ; polygon(ref.line,col="white")
f.image(varName='nonStatKriSD',polyName='myPoly',col=myPalette, zlim=c(0,5)) ; polygon(ref.line,col="white")
f.image(varName='nonStatKriCV',polyName='myPoly',col=myPalette,zlim=c(0,0.4)); polygon(ref.line,col="white")
dev.print(device = png, file = "Res/CVNonStat.png",width=450,height=450)

f.image(varName='nonStatSimu1',polyName='myPoly',zlim=c(0,300))
f.image(varName='nonStatSimu2',polyName='myPoly',zlim=c(0,300))
f.image(varName='nonStatSimu3',polyName='myPoly',zlim=c(0,300))
f.image(varName='nonStatSimu4',polyName='myPoly',zlim=c(0,300))

image(unique(myGrid['x1']), unique(myGrid['x2']),
      matrix(ifelse(myGrid['Polygon']==1 & abs(myGrid['nonStatKriCV']) < 0.2,myGrid['nonStatKri'],NA),nx),
      zlim=c(0,200),
      asp=1,col=myPalette,
      xlab="Longitude",ylab="Latitude",las=1,
      xlim=range(myGrid['x1'][myGrid['Polygon']==1]))
lines(ref.line)

tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("nonStatSimu",1)],NA),nx))
rowRank <- 125
f.image(varName='nonStatKri',zlim=c(0,100),las=1) ; polygon(ref.line,col="white")
abline(v=tempRes$x[rowRank])
plot(tempRes$y, -tempRes$z[rowRank,],type="l",las=1,ylab="Depth",xlab="Latitude",
     ylim=c(-70,0),xlim=c(-3.92,-3.77))
for(i in 2:10){
  tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("nonStatSimu",i)],NA),nx))
lines(tempRes$y,-tempRes$z[rowRank,],col=i)
}

dev.print(device = png, file = "Res/simuCondNonStat.png",width=600,height=400)

                
