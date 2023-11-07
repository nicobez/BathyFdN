########################################
# SPDE Stationary kriging (KO? KS?) with the ML parameters myOptimStatParamStd
########################################

sel <- myGrid['delatedPolygon']
selRank <- (1:length(sel))[sel==1]
temp <- rep(NA,length(sel))

ndat <- length(myDb['depthStd'])

myModelStat = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=myOptimStatParamStd[c(2,3)],sill=1)

# building D and invD
nugget = myOptimStatParamStd[1]
DStat <- Diagonal(ndat,nugget)
invDStat <- Diagonal(ndat,1/nugget)

# building Q
QStat = cs_toTL(PrecisionOpCs(myMesh,myModelStat)$getQ()) #Precision matrix for sill=1

# Solving 
M1 <- QStat + t(A) %*% invDStat %*% A
M2 <- t(A) %*% invDStat %*% myDb["depthStd"] 
zKri <-  Matrix:::solve(M1,M2)

# Affecting the result depth = depthStd*sd+m to the grid
temp[selRank] <-zKri
myGrid['statKri'] <- driftSd(myGrid['distAcross'])*temp + driftM(myGrid['distAcross'])

###########################################
# Kriging estimation variance Using NonCond simulations
###########################################
cholQStat = chol(QStat)

nSimu <- 400 # 40 secondes for nSimu = 200
simuCond <- array(NA,dim=c(myGrid$getActiveSampleNumber(),nSimu))
for(i in 1:nSimu){
  simuNonCondVertices = solve(cholQStat,rnorm(NROW(QStat),0,1))[,1]
  simuNonCondData <- A %*% simuNonCondVertices + rnorm(NROW(A),0,sqrt(nugget))
  M2 <- t(A) %*% invDStat %*% simuNonCondData 
  nonCondKri <- Matrix:::solve(M1,M2)
  simuCondVertices <- zKri + simuNonCondVertices - nonCondKri
  simuCond[,i] <- as.vector(simuCondVertices)
}

# depth = depthStd*sd + m
simuCond <- driftSd(myGrid['distAcross'][selRank])*simuCond + driftM(myGrid['distAcross'][selRank])
# truncation of depth to the observed extremum
simuCond[simuCond < min(myDb['depth'],na.rm=T)] <- min(myDb['depth'],na.rm=T)

# Affecting the results to the grid
sdKriDepth <- apply(simuCond,1,"sd")
temp[selRank] <- sdKriDepth
sdKriDepth <- temp

for(i in 1:10){
  temp[selRank] <- simuCond[,i]
  assign(paste0("statSimu",i),temp)
}

myGrid['statKriSD'] <- sdKriDepth
myGrid['statKriCV'] <- sdKriDepth/myGrid['statKri']
for(i in 1:10) myGrid[paste0("statSimu",i)] <- get(paste0("statSimu",i))

f.image(varName='statKri',polyName='myPoly',zlim=c(0,200)) ; polygon(ref.line,col="white")
points(myDb['longitude'],myDb['latitude'],cex=0.1)
lines(polyTotal)

f.image(varName='statKriSD',polyName='myPoly',col=myPalette, zlim=c(0,5))
f.image(varName='statKriCV',polyName='myPoly',col=myPalette,zlim=c(0,0.2)) ; polygon(ref.line,col="white") ; lines(polyTotal)
f.image(varName='statSimu1',polyName='myPoly',zlim=c(0,300))
f.image(varName='statSimu2',polyName='myPoly',zlim=c(0,300))
f.image(varName='statSimu3',polyName='myPoly',zlim=c(0,300))
f.image(varName='statSimu4',polyName='myPoly',zlim=c(0,300))

image(unique(myGrid['x1']), unique(myGrid['x2']),
      matrix(ifelse(myGrid['Polygon']==1 & abs(myGrid['statKriCV']) < 0.2,myGrid['statKri'],NA),nx),
      zlim=c(0,200),
      asp=1,col=myPalette,
      xlab="Longitude",ylab="Latitude",las=1,
      xlim=range(myGrid['x1'][myGrid['Polygon']==1]))
lines(ref.line)


tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("statSimu",1)],NA),nx))
rowRank <- 125
f.image(varName='statKri',zlim=c(0,100))
abline(v=tempRes$x[rowRank])
plot(-tempRes$z[rowRank,],ylim=c(-100,0),type="l")
for(i in 2:10){
  tempRes <- list(x=unique(myGrid['x1']), 
                y=unique(myGrid['x2']),
                z=matrix(ifelse(myGrid['Polygon']==1,myGrid[paste0("statSimu",i)],NA),nx))
lines(-tempRes$z[rowRank,],col=i)
}
                

