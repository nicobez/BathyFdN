# Empirical conditional expectations of order 1 & 2 ----------------------
# Computing the drift function with the HD data
x <- myDb['distAcross']
z <- myDb['depth']

#x.breaks <- seq(0,max(x),0.001)
distAcrossBreaks <- c(seq(0,1.2,length=50),seq(1.2,1.5,length=5))[-50]#max(x,na.rm=T),length=10))[-100]
distAcrossMids <- (distAcrossBreaks[-1]+distAcrossBreaks[-length(distAcrossBreaks)])/2

tempMean <- tapply(z,cut(x,breaks=distAcrossBreaks),FUN="mean")
tempSd <- tapply(z,cut(x,breaks=distAcrossBreaks),FUN="sd")
driftM <- approxfun(x=distAcrossMids,y=tempMean,rule=2)
driftSd <- approxfun(x=distAcrossMids,y=tempSd,rule=2)

# Cleaning
rm(x,z,tempMean,tempSd)

# Compute centered and stadardized depth -------------------------------------------------------
depthC <- (myDb['depth']-driftM(myDb['distAcross']))
depthStd <- (myDb['depth']-driftM(myDb['distAcross']))/driftSd(myDb['distAcross'])
myDb$addColumns(depthC,'depthC')
myDb$addColumns(depthStd,'depthStd')

# residuals <- (depthHighDef['depth']-driftM(depthHighDef['distAcross']))/driftSd(depthHighDef['distAcross'])
# depthHighDef$addColumns(residuals,'depthResiduals',ELoc_Z(),0)

rm(depthC,depthStd)
