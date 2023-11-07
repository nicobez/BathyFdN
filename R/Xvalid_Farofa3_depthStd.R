
### db containing the data points IN the criteria to be excluded, 'in' the polygon
xValidIn <- myDb$clone()
xValidIn$addSelection(xValidIn['farofa']==3)
xValidIn <- Db_createReduce(xValidIn)
xValidIn$deleteColumns('NewSel')

### db containing the data points OUT the criteria to be excluded, 'out' the polygon
xValidOut <- myDb$clone()
xValidOut$addSelection(xValidOut['farofa']<3)
xValidOut <- Db_createReduce(xValidOut)
xValidOut$deleteColumns('NewSel')

### Nota using *$deleteSamples is much too longer

# krigings without the crossValidation points
# QStat, QNoStat and myMesh are unchanged
AOut = cs_toTL(ProjMatrix_create(xValidOut,myMesh)$getAproj())

nOut <- xValidOut$getSampleNumber()

# Stationary kriging
DStat <- Diagonal(nOut,nugget)
invDStat <- Diagonal(nOut,1/nugget)
M1 <- QStat + t(AOut) %*% invDStat %*% AOut
M2 <- t(AOut) %*% invDStat %*% xValidOut["depthStd"] 
zKriStat <-  Matrix:::solve(M1,M2)
#zKriStat <- driftSd(myGrid['distAcross'][selRank])*zKriStat + driftM(myGrid['distAcross'][selRank])

# NonStationary kriging
DNoStat <- Diagonal(nOut,xValidOut['nugget'])
invDNoStat <- Diagonal(nOut,1/xValidOut['nugget'])
M1 <- QNoStat + t(AOut) %*% invDNoStat %*% AOut
M2 <- t(AOut) %*% invDNoStat %*% xValidOut["depthStd"] 
zKriNoStat <-  Matrix:::solve(M1,M2)
#zKriNoStat <- driftSd(myGrid['distAcross'][selRank])*zKriNoStat + driftM(myGrid['distAcross'][selRank])

# # affecting the result to the grid
# temp[selRank] <- zKriStat
# myGrid['statkriXValid'] <- temp
# 
# f.image(varName='kriXValid',zlim=c(0,200))

# kriging at crossValidation points :
# projection of the vertices' values towards the crossValidation points
AIn = cs_toTL(ProjMatrix_create(xValidIn,myMesh)$getAproj())

xValidStat <- as.numeric(AIn %*% zKriStat)
xValidNoStat <- as.numeric(AIn %*% zKriNoStat) 

#z <- as.numeric(xValidIn['depth'])
z <- as.numeric(xValidIn['depthStd'])

xValidIn$addColumns(xValidNoStat,'xValidNoStat')
xValidIn$addColumns(xValidStat,'xValidStat')
xValidIn$addColumns(xValidStat-z,'xValidErrorStat')
xValidIn$addColumns(xValidNoStat-z,'xValidErrorNoStat')
xValidIn$addColumns(xValidStat-xValidNoStat,'xValidDiff')

ggplot() + plot(xValidIn,name_color="xValidDiff")


plot(z,xValidNoStat,ylim=range(c(xValidNoStat,xValidStat)),pch=20,col=rgb(0,0,0,0.5),ylab="z.xValid")
abline(0,1,col=3,lwd=3)
points(z,xValidStat,pch=20,col=rgb(1,0,0,0.5))
segments(z,xValidNoStat,z,xValidStat,lty=3)

plot(z,xValidNoStat,pch=20,col=rgb(0,0,0,0.5),ylab="z.xValid",xlim=c(0,80),ylim=c(0,80))
abline(0,1,col=3,lwd=3)
points(z,xValidStat,pch=20,col=rgb(1,0,0,0.5))
segments(z,xValidNoStat,z,xValidStat,lty=3)

plot(z,xValidNoStat-z,ylim=range(c(xValidNoStat-z,xValidStat-z)),pch=20,col=rgb(0,0,0,0.5))
abline(h=0,col=3,lwd=3)
points(z,xValidStat-z,pch=20,col=rgb(1,0,0,0.5))
segments(z,xValidNoStat-z,z,xValidStat-z,lty=3)

plot(z[xValidIn['distAcross']<1.1],(xValidNoStat-xValidStat)[xValidIn['distAcross']<1.1])
abline(h=0,col=2,lwd=2);abline(0,1,col=2) ; abline(0,-1,col=2)

plot(xValidIn['distAcross'],xValidNoStat-z,xlim=c(0,0.9),ylim=c(-2,2))
points(xValidIn['distAcross'],xValidStat-z,col=rgb(1,0,0,0.5))


hist(xValidNoStat-xValidStat,nclass=50)

plot(ecdf(xValidNoStat))
plot(ecdf(xValidStat),add=T,col=2)
abline(v=0,col=3)

plot(ecdf((xValidNoStat-xValidStat)))

sel <- (xValidNoStat-xValidStat) > 60
xValidIn$addSelection(sel)


