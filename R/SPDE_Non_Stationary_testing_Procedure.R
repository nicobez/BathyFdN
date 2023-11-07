f.image <- function(grid=myGrid,varName,col=rev(myPalette),polyName=NULL,...){
  nx <- grid$getNX(0)
  image(unique(grid['x1']), unique(grid['x2']),
        if(is.null(polyName)) matrix(grid[varName],nx) 
        else matrix(ifelse(grid['Polygon']==1,grid[varName],NA),nx),
        col=col,
        xlab="Longitude",ylab="Latitude",las=1,...)
}


###########
# Testing procedure

# Building the grid
testGrid <- DbGrid_create(x0=c(-2,-2),
                        dx=c(0.1,0.1)/4,
                        nx=c(40,40)*4)
testGrid['dir'] <- atan(testGrid['x2']/testGrid['x1'])/pi*180
testGrid['dist'] <- sqrt(testGrid['x1']^2 + testGrid['x2']^2)
testGrid['r1'] <- ifelse(testGrid['dist'] < 1,0.2 + (1-0.2)*testGrid['dist'],1)
testGrid['r2'] <- ifelse(testGrid['dist'] < 1,1 + (3-1)*testGrid['dist'],3)
testGrid$setLocator("dir",ELoc_NOSTAT(),0)
testGrid$setLocator("r1",ELoc_NOSTAT(),1)
testGrid$setLocator("r2",ELoc_NOSTAT(),2)

# buidling the Db with non stationary nugget
testDb <- Db_create()
testDb$addColumns(c(-0.5,1),'longitude',ELoc_X(),0)
testDb$addColumns(c(-0.5,1),'latitude',ELoc_X(),1)
testDb$addColumns(c(3,3),'z',ELoc_Z(),0)
#testDb$addColumns(c(0.1,0.5),'nugget') 
testDb$addColumns(c(0.1,0.1),'nugget') 

# building the non stationary model with sill=1
nostat = NoStatArray(c("A1","R1","R2"),testGrid)
testModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=c(1,1),sill=1)
Model_addNoStat(testModel,nostat)

# building the key SPDE matrices
testMesh <- MeshETurbo(testGrid)
testA <- cs_toTL(ProjMatrix_create(testDb,testMesh)$getAproj())
testQ <- cs_toTL(PrecisionOpCs(testMesh,testModel)$getQ())
testD <- Diagonal(2,testDb['nugget'])
testinvD <- Diagonal(2,1/testDb['nugget'])

# Kriging
testM1 <- testQ + t(testA) %*% testinvD %*% testA
testM2 <- t(testA) %*% testinvD %*% testDb["z"] 
testKri <-  Matrix:::solve(testM1,testM2)
testGrid['kriging'] <- as.numeric(testKri)

f.image(grid=testGrid,varName='kriging',col=myPalette)

# Conditional Simulations
testcholQ = chol(testQ)

nSimu <- 400
testsimuCond <- array(NA,dim=c(testGrid$getActiveSampleNumber(),nSimu))
for(i in 1:nSimu){
  simuNonCondVertices = solve(testcholQ,rnorm(NROW(testQ),0,1))[,1]
  simuNonCondData <- testA %*% simuNonCondVertices + rnorm(NROW(testA),0,sqrt(testDb['nugget']))
  testM2 <- t(testA) %*% testinvD %*% simuNonCondData 
  nonCondKri <- Matrix:::solve(testM1,testM2)
  simuCondVertices <- testKri + simuNonCondVertices - nonCondKri
  testsimuCond[,i] <- as.numeric(simuCondVertices)
}

# Affecting outputs to the grid
meanSimu <- apply(testsimuCond,1,"mean")
sdSimu <- apply(testsimuCond,1,"sd")
testGrid['nonStatKriSD'] <- sdSimu
testGrid['nonStatKriCV'] <- sdSimu/testGrid['kriging']
for(i in 1:10) testGrid[paste0("Simu",i)] <- testsimuCond[,i]

plot(testKri,meanSimu);abline(0,1,col=2,lwd=2) # to check that mean of N cond simu = kriging
f.image(grid=testGrid,varName='kriging',col=myPalette)
f.image(grid=testGrid,varName='nonStatKriSD',col=myPalette)
f.image(grid=testGrid,varName='nonStatKriCV',col=myPalette,zlim=c(0,0.5))
f.image(grid=testGrid,varName='Simu1',col=myPalette)
f.image(grid=testGrid,varName='Simu5',col=myPalette)



# buidling a larger Db
set.seed(12345)
testDb <- Db_create()
testDb$addColumns(sort(runif(25,-2,2)),'longitude',ELoc_X(),0)
testDb$addColumns(sort(runif(25,-2,2)),'latitude',ELoc_X(),1)
testDb$addColumns(1:25,'z',ELoc_Z(),0)

# buidling a coarser Db
testGrid <- DbGrid_create(x0=c(-2,-2),
                          dx=c(0.25,0.25),
                          nx=c(20,20))

toto <- dbStatisticsPerCell(testDb,testGrid,oper=EStatOption(1),name1='z',name2='z')
testGrid$addColumns(toto,'z')

ggplot() + plot(testDb,name_color="z")

ggplot() + plot(testGrid,name_raster="z")



