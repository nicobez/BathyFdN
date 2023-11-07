#myDb$setLocator('distAcross',ELoc_X(),0)
#myDb$setLocator('distAlong',ELoc_X(),1)

if(compute){
  myDb$setLocator('depthStd',ELoc_Z(),0)
  vgY <- list()
  k <- 1
  for(i in seq(0,1.2,0.1)){ 
    cat("\n",k,"\n")
    tempDb <- myDb$clone()
    selDist <- tempDb['distAcross'] > i & tempDb['distAcross'] < (i+0.1)
    tempDb$addSelection(selDist)
    
    tempVarioParamOmni = VarioParam()
    tempDir = DirParam_create(npas=30,dpas=pixelSize)
    tempVarioParamOmni$addDir(tempDir)
    
    tempVarioOmni = Vario(tempVarioParamOmni,tempDb)
    tempVarioOmni$compute(ECalcVario_VARIOGRAM())
    
    resg <- tempVarioOmni$getAllGg(0)
    resh <- tempVarioOmni$getAllHh(0)
    resw <- tempVarioOmni$getAllSw(0)
    
    vgY[[k]] <- list(h=resh,g=resg,w=resw)
    #assign(paste0("vgY",k),res)
    k <- k+1
  }
  
  myDb$setLocator('depth',ELoc_Z(),0)
  k <- 1
  vgZ <- list()
  for(i in seq(0,1.2,0.1)){ 
    cat("\n",k,"\n")
    tempDb <- myDb$clone()
    selDist <- tempDb['distAcross'] > i & tempDb['distAcross'] < (i+0.1)
    tempDb$addSelection(selDist)
    
    tempVarioParamOmni = VarioParam()
    tempDir = DirParam_create(npas=30,dpas=pixelSize)
    tempVarioParamOmni$addDir(tempDir)
    
    tempVarioOmni = Vario(tempVarioParamOmni,tempDb)
    tempVarioOmni$compute(ECalcVario_VARIOGRAM())
    
    resg <- tempVarioOmni$getAllGg(0)
    resh <- tempVarioOmni$getAllHh(0)
    resw <- tempVarioOmni$getAllSw(0)
    
    vgZ[[k]] <- list(h=resh,g=resg,w=resw)
    #assign(paste0("vgZ",k),res)
    k <- k+1
  }
  
  dump('vgY',file="Data/vgY")
  dump('vgZ',file="Data/vgZ")
}

source('Data/vgY')
source('Data/vgZ')

#########

vg <- vgY #vgY

plot(vg[[1]]$h,vg[[1]]$g, type="b",lwd=2,xlim=c(0,0.03),xaxs="i",ylim=c(0,2),yaxs="i") ; abline(h=1)
for(i in 2:12) {
  lines(vg[[i]]$h,vg[[i]]$g,col=i,type="b",lwd=2)
  browser()
}


plot(seq(0,1.1,0.1)+0.05,log10(sapply(1:12, function(i) mean(vg[[i]]$g[3:5])/mean(vg[[i]]$h[3:5]))),type='b')
plot(seq(0,0.8,0.1)+0.05,sapply(1:9, function(i) mean(vg[[i]]$g[3:5])/mean(vg[[i]]$h[3:5])),type='b')

symbols(seq(0,1.1,0.1)+0.05,
        sapply(1:12, function(i) mean(vg[[i]]$h[1:5])/mean(vg[[i]]$g[1:5])),
        sqrt(sapply(1:12, function(i) mean(vg[[i]]$w[1:5]))),
        inches=0.5)
lines(seq(0,1.1,0.1)+0.05,sapply(1:12, function(i) mean(vg[[i]]$h[1:5])/mean(vg[[i]]$g[1:5])))

setAnisoRange2 <- approxfun(x=seq(0,1.1,0.1)+0.05,y=sapply(1:12, function(i) mean(vg[[i]]$h[1:5])/mean(vg[[i]]$g[1:5])),rule=2)

lines(seq(0,1.3,0.05),setAnisoRange2(seq(0,1.3,0.05)),col=2)

symbols(seq(0,1.2,0.1)+0.05,sapply(1:13, function(i) get(paste0("vg",i))$g[3]),circles=sqrt(sapply(1:13, function(i) get(paste0("vg",i))$w[3])))
lines(seq(0,1.2,0.1)+0.05,sapply(1:13, function(i) get(paste0("vg",i))$g[3]))


plot(seq(0,1.2,0.1)+0.05,sapply(1:13, function(i) mean(get(paste0("vg",i))$g[5:7])),type='b')





myvarioParam = VarioParam()
mydirs = DirParam_createMultiple(ndir=2, npas=100, dpas=pixelSize)
myvarioParam$addMultiDirs(mydirs)
myvario = Vario(myvarioParam,myDb)
myvario$compute(ECalcVario_VARIOGRAM())


if (verbose)
  myvario$display()

if (graphics)
  ggplot() + plot.varmod(myvario) + 
  plot.decoration(title="Multi-Directional Variogram of Pb")
