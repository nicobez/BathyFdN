
set.seed(1234)
selxValidRank <- sample(0:ndat, 5000) # use sample(0:(ndat-1), 5000) if the selection is made on the ranks

selxValidOut <- rep(T,ndat)
selxValidOut[selxValidRank] <- F
selxValidIn <- !selxValidOut

xValidIn <- myDb$clone()
xValidOut <- myDb$clone()

xValidIn$addSelection(selxValidIn)
xValidIn = Db_createReduce(xValidIn)
xValidIn$deleteColumns('NewSel')

xValidOut$addSelection(selxValidOut)
xValidOut = Db_createReduce(xValidOut)
xValidOut$deleteColumns('NewSel')

### Nota using *$deleteSamples is much too longer

# krigings without the crossValidation points
# QStat and QNoStat are unchanged
AprojOut = cs_toTL(ProjMatrix_create(xValidOut,myMesh)$getAproj())

nOut <- xValidOut$getSampleNumber()

# Stationary kriging
DStat <- Diagonal(nOut,nugget)
invDStat <- Diagonal(nOut,1/nugget)
M1 <- QStat + t(AprojOut) %*% invDStat %*% AprojOut
M2 <- t(AprojOut) %*% invDStat %*% xValidOut["depthStd"] 
zKriStat <-  Matrix:::solve(M1,M2)
zKriStat <- driftSd(myGrid['distAcross'][selRank])*zKriStat + driftM(myGrid['distAcross'][selRank])

# NonStationary kriging
DNoStat <- Diagonal(nOut,xValidOut['nugget'])
invDNoStat <- Diagonal(nOut,1/xValidOut['nugget'])
M1 <- QNoStat + t(AprojOut) %*% invDNoStat %*% AprojOut
M2 <- t(AprojOut) %*% invDNoStat %*% xValidOut["depthStd"] 
zKriNoStat <-  Matrix:::solve(M1,M2)
zKriNoStat <- driftSd(myGrid['distAcross'][selRank])*zKriNoStat + driftM(myGrid['distAcross'][selRank])

# # affecting the result to the grid
# temp[selRank] <- zKriStat
# myGrid['statkriXValid'] <- temp
# 
# f.image(varName='kriXValid',zlim=c(0,200))

# kriging at crossValidation points :
# projection of the vertices' values towards the crossValidation points
AprojIn = cs_toTL(ProjMatrix_create(xValidIn,myMesh)$getAproj())

xValidStat <- as.numeric(AprojIn %*% zKriStat)
xValidNoStat <- as.numeric(AprojIn %*% zKriNoStat) 

z <- as.numeric(xValidIn['depth'])

plot(z,xValidNoStat,ylim=range(c(xValidNoStat,xValidStat)))
abline(0,1,col=2,lwd=2)
points(z,xValidStat,col=2)
segments(z,xValidNoStat,z,xValidStat,lty=2)

plot(xValidNoStat,xValidStat)
abline(0,1,col=2,lwd=2)

hist(xValidNoStat-xValidStat)
