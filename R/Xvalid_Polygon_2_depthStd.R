
xValiPolygonDf <- read.table('Data/xvalid.polygon.4.ascii',skip=1)
temp <- Db_create()
temp$addColumns(tab=xValiPolygonDf[,1],"longitude",ELoc_X(),0)
temp$addColumns(tab=xValiPolygonDf[,2],"latitude",ELoc_X(),1)
xValiPolygon <- Polygons_createFromDb(temp)

xValidIn <- myDb$clone()
db_polygon(xValidIn, xValiPolygon)
temp <- !xValidIn['Polygon']
xValidIn <- Db_createReduce(xValidIn)
xValidIn$deleteColumns('Polygon')

xValidOut <- myDb$clone()
xValidOut$addColumns(temp,'Polygon',ELoc_SEL(),0)
xValidOut <- Db_createReduce(xValidOut)
xValidOut$deleteColumns('Polygon')
               
f.image(varName='nonStatKriCV',col=myPalette,zlim=c(0,0.2))
lines(xValiPolygonDf,lwd=2)

# selCrossValidRank <- sample(0:ndat, 5000) # use sample(0:(ndat-1), 5000) if the selection is made on the ranks
# 
# selCrossValidIn <- rep(T,ndat)
# selCrossValidIn[selCrossValidRank] <- F
# selCrossValidOut <- !selCrossValidIn
# 
# myDbCrossValidIn <- myDb$clone()
# myDbCrossValidOut <- myDb$clone()
# 
# myDbCrossValidIn$addSelection(selCrossValidIn)
# myDbCrossValidIn = Db_createReduce(myDbCrossValidIn)
# myDbCrossValidIn$deleteColumns('NewSel')
# 
# myDbCrossValidOut$addSelection(selCrossValidOut)
# myDbCrossValidOut = Db_createReduce(myDbCrossValidOut)
# myDbCrossValidOut$deleteColumns('NewSel')

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

plot(z,xValidNoStat,ylim=range(c(xValidNoStat,xValidStat)),pch=20,col=rgb(0,0,0,0.5),ylab="z.xValid")
abline(0,1,col=3,lwd=3)
points(z,xValidStat,pch=20,col=rgb(1,0,0,0.5))
segments(z,xValidNoStat,z,xValidStat,lty=3)


plot(z,xValidNoStat-z,ylim=range(c(xValidNoStat-z,xValidStat-z)),pch=20,col=rgb(0,0,0,0.5))
abline(h=0,col=3,lwd=3)
points(z,xValidStat-z,pch=20,col=rgb(1,0,0,0.5))
segments(z,xValidNoStat-z,z,xValidStat-z,lty=3)

plot(xValidNoStat,xValidStat)
abline(0,1,col=2,lwd=2)

hist(xValidNoStat-xValidStat,nclass=50)
