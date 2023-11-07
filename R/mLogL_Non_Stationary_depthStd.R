
####################### ML estimation of the non stationary parameters
# 
# Estimation of the non stationary range parameter values and the nugget by max Likelihood
# Instead of doing an optimization, we optimize the likelihood on progressively refined experimental designs.
# The result could serve as relevant starting points for an eventual optim() run after.
# (not done .... the "price" of doing it is not obvious) 
# 
##############################################

myModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=c(1,1),sill=1)

#############################
# Coarse experimental design (~2 hours)
#############################
if(compute){
  nuggetStart <- c(0.01,0.05,0.1)
  nuggetEnd <- c(0.01,0.05,0.1)
  r1Start <- c(0.01,0.025,0.05)
  r1End <- c(0.01,0.025,0.05)
  r2Start <- c(0.01,0.025,0.05)
  r2End <- c(0.01,0.025,0.05)
  planExpStd1 <- expand.grid(nuggetStart,nuggetEnd,r1Start,r1End,r2Start,r2End)
  planExpStd1[,7] <- rep(NA,NROW(planExpStd1))
  names(planExpStd1) <- c("nuggetStart","nuggetEnd","r1Start","r1End","r2Start","r2End","mLogL")
  for(i in 1:NROW(planExpStd1)){
    cat("\n", i,  "/", NROW(planExpStd1),"\n")
    temp <- mLogLDepthStd(optimParam=as.numeric(planExpStd1[i,1:6]),dat=myDb["depthStd"],
                          model=myModel,Aproj=Aproj,mesh=myMesh,
                          mean=0,fixParam=1,NoStat = T,grid=myGrid,db=myDb)
    planExpStd1[i,7] <- temp
  }
  planExpStd1 <- planExpStd1[order(planExpStd1[,7]),]
  save("planExpStd1",file='Res/planExpStd1')
} else{
  load('Res/planExpStd1')
}

plot(planExpStd1[,7],type="b")
plot(planExpStd1[,7][1:20],type="b")

plot(planExpStd1[,3],planExpStd1[,7])

# Comments:
# the minimum LogL correspond all to nugStart = nugEnd = 0.01
# These two parameters can be fixed to these values.
# ==> nugStart = 0.01, nugEnd = 0.01
# The ranges in the across direction decreases 
# The ranges in the along direction alternate between 0.01 and 0.025

###############################
# Refined experimental design around the optimal point of the former experimental design for the ranges in the along direction
##############################
if(compute){
  nuggetStart <- planExpStd1[1,1]
  nuggetEnd <- planExpStd1[1,2]
  r1Start <- planExpStd1[1,3] + c(-0.005,0,0.005)
  r1End <- planExpStd1[1,4] + c(-0.005,0,0.005)
  r2Start <- planExpStd1[1,5] + c(-0.005,0,0.005) #c(0.005,0.01,0.015,0.02,0.025,0.03)
  r2End <- planExpStd1[1,6]+ c(-0.005,0,0.005) #c(0.005,0.01,0.015,0.02,0.025,0.03)
  planExpStd2 <- expand.grid(nuggetStart,nuggetEnd,r1Start,r1End,r2Start,r2End)
  planExpStd2[,7] <- rep(NA,NROW(planExpStd2))
  names(planExpStd2) <- c("nuggetStart","nuggetEnd","r1Start","r1End","r2Start","r2End","mLogL")
    for(i in 1:NROW(planExpStd2)){
    cat("\n", i, "/", NROW(planExpStd2), "\n")
    temp <- mLogLDepthStd(optimParam=as.numeric(planExpStd2[i,1:6]),dat=myDb["depthStd"],
                  model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                  NoStat = T,grid=myGrid,db=myDb)
    planExpStd2[i,7] <- temp
    }
  planExpStd2 <- planExpStd2[order(planExpStd2[,7]),]
  save("planExpStd2",file='Res/planExpStd2')
} else{
  load('Res/planExpStd2')
}

plot(planExpStd2[,7],type="b")

# This 2d experimental designa confirms that the ranges decrease in both directions when going off-shore
# r1Start is optimum for a medium value ==> kept
# the other range values converge towards the borders of the intervals ==> should be re analysed

###############################
# Refined experimental design ... continue
##############################
if(compute){
  nuggetStart <- planExpStd1[1,1]
  nuggetEnd <- planExpStd1[1,2]
  r1Start <- planExpStd2[1,3] 
  r1End <- planExpStd2[1,4] + c(-0.0025,0,0.0025)
  r2Start <- planExpStd2[1,5] + c(-0.0025,0,0.0025) #c(0.005,0.01,0.015,0.02,0.025,0.03)
  r2End <- planExpStd2[1,6]+ c(-0.0025,0,0.0025) #c(0.005,0.01,0.015,0.02,0.025,0.03)
  planExpStd3 <- expand.grid(nuggetStart,nuggetEnd,r1Start,r1End,r2Start,r2End)
  planExpStd3[,7] <- rep(NA,NROW(planExpStd3))
  names(planExpStd3) <- c("nuggetStart","nuggetEnd","r1Start","r1End","r2Start","r2End","mLogL")
  for(i in 1:NROW(planExpStd3)){
    cat("\n", i, "/", NROW(planExpStd3), "\n")
    temp <- mLogLDepthStd(optimParam=as.numeric(planExpStd3[i,1:6]),dat=myDb["depthStd"],
                          model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                          NoStat = T,grid=myGrid,db=myDb)
    planExpStd3[i,7] <- temp
  }
  planExpStd3 <- planExpStd3[order(planExpStd3[,7]),]
  save("planExpStd3",file='Res/planExpStd3')
} else{
  load('Res/planExpStd3')
}

plot(planExpStd3[,7],type="b")

###############################
# Refined experimental design ... continue
##############################
if(compute){
  nuggetStart <- planExpStd3[1,1]
  nuggetEnd <- planExpStd3[1,2]
  r1Start <- planExpStd3[1,3] 
  r1End <- planExpStd3[1,4] + c(-0.0025,0,0.0025)
  r2Start <- planExpStd3[1,5] + c(-0.0025,0,0.0025) #c(0.005,0.01,0.015,0.02,0.025,0.03)
  r2End <- planExpStd3[1,6]
  planExpStd4 <- expand.grid(nuggetStart,nuggetEnd,r1Start,r1End,r2Start,r2End)
  planExpStd4[,7] <- rep(NA,NROW(planExpStd4))
  names(planExpStd4) <- c("nuggetStart","nuggetEnd","r1Start","r1End","r2Start","r2End","mLogL")
  for(i in 1:NROW(planExpStd4)){
    cat("\n", i, "/", NROW(planExpStd4), "\n")
    temp <- mLogLDepthStd(optimParam=as.numeric(planExpStd4[i,1:6]),dat=myDb["depthStd"],
                          model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                          NoStat = T,grid=myGrid,db=myDb)
    planExpStd4[i,7] <- temp
  }
  planExpStd4 <- planExpStd4[order(planExpStd4[,7]),]
  save("planExpStd4",file='Res/planExpStd4')
} else{
  load('Res/planExpStd4')
}

plot(planExpStd4[,7],type="b")


#######################################
#######################################
#######################################
myOptimNoStatParamStd = as.numeric(planExpStd4[1,1:6])



