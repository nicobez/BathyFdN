
####################### ML estimation of the parameters
# Estimation of the non stationary range parameter values and the proportion of nugget by max Likelihood
# Before doing an optimization we search for relevant starting points based on 
# a coarse experimental design  
# and a refined one around the minimum obtain with the coarse design
# 
# An optim() run is possible after (not done .... the "price" of doing it is not obvious)
##############################################


myModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=c(1,1),sill=1)

# First trial and timer
myOptimParamInit <- c(0.05,0.03,0.01,0.02,0.015)
start <- Sys.time()
mLogL0 = abs(mLogL(myOptimParamInit,myDb["depthC"],myModel,Aproj,myMesh,NoStat = T,grid=myGrid,db=myDb,mean=0))
print( Sys.time() - start )
# ~ 8 seconds per logL

#############################
# Coarse experimental design
#############################
propNuggetPlanExp <- c(0.01,0.05,0.1)
r1StartPlanExp <- c(0.01,0.025,0.05)
r1EndPlanExp <- c(0.025,0.05,0.75)
r2StartPlanExp <- c(0.01,0.025,0.05)
r2EndPlanExp <- c(0.025,0.05,0.75)
planExp <- expand.grid(propNuggetPlanExp,r1StartPlanExp,r1EndPlanExp,r2StartPlanExp,r2EndPlanExp)
names(planExp) <- c("propNugget","r1Start","r1End","r2Start","r2End")

if(compute){
  planExpLogl <- rep(NA,NROW(planExp))
  for(i in 1:NROW(planExp)){
    cat(i, "\n")
    tempParam <- as.numeric(planExp[i,]) # c(0.02,as.numeric(planExp[i,])) # tempParam <- c(0.02,as.numeric(planExp[1,]))
    temp <- mLogL(optimParam=tempParam,dat=myDb["depthC"],model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                  NoStat = T,grid=myGrid,db=myDb)
    print(temp)
    planExpLogl[i] <- temp
  }
  save("planExpLogl",file='Res/planExpLogl')
} else{
  load('Res/planExpLogl')
}

plot(sort(planExpLogl),type="b")
plot(sort(planExpLogl)[1:20],type="b")
coarseRes <- planExp[order(planExpLogl),]
coarseOptim <- as.numeric(coarseRes[1,])
plot(coarseRes[,1],planExpLogl[order(planExpLogl)])


# Comments:
# the first sets of parameters correspond all to propNugget = 0.1, then 0.05 and 0.01
# the first four sets of parameters argue for 
#   an increasing range in the across direction from coast to offshore r1Start < r1End
#   increasing range in the along direction 

###############################
# Refined experimental design around the optimal point of the coarse one
# except for the propNugget that seems already well estimated
##############################
r1StartPlanExp2 <- coarseOptim[2] + c(-0.01,0,0.01)
r1EndPlanExp2 <- coarseOptim[3] + c(-0.01,0,0.01)
r2StartPlanExp2 <- coarseOptim[4] + c(-0.01,0,0.01)
r2EndPlanExp2 <- coarseOptim[5] + c(-0.01,0,0.01)
planExp2 <- expand.grid(r1StartPlanExp2,r1EndPlanExp2,r2StartPlanExp2,r2EndPlanExp2)
names(planExp2) <- c("r1Start","r1End","r2Start","r2End")

if(compute){
  planExpLogl2 <- rep(NA,NROW(planExp2))
  for(i in 1:NROW(planExp2)){
    cat(i, "\n")
    tempParam <- c(coarseOptim[1],as.numeric(planExp2[i,])) # tempParam <- c(0.02,as.numeric(planExp[1,]))
    temp <- mLogL(optimParam=tempParam,dat=myDb["depthC"],model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                  NoStat = T,grid=myGrid,db=myDb)
    print(temp)
    planExpLogl2[i] <- temp
  }
  save("planExpLogl2",file='Res/planExpLogl2')
} else{
  load('Res/planExpLogl2')
}

plot(sort(planExpLogl2),type="b")
plot(sort(planExpLogl2)[1:20],type="b")
refinedRes <- planExp2[order(planExpLogl2),]
refinedOptim <- as.numeric(refinedRes[1,])
plot(refinedRes[,4],planExpLogl2[order(planExpLogl2)])

# the r1start and r2start are unchanged
# instead r1End and r2End are still equal to the upper right limit

###############################
# Refined experimental design around the optimal point for r1End and r2End
##############################
r1EndPlanExp3 <- refinedOptim[2] + c(0,0.01,0.02)
r2EndPlanExp3 <- refinedOptim[4] + c(0,0.01,0.02)
planExp3 <- expand.grid(r1EndPlanExp3,r2EndPlanExp3)
names(planExp3) <- c("r1End","r2End")

if(compute){
  planExpLogl3 <- rep(NA,NROW(planExp3))
  for(i in 1:NROW(planExp3)){
    cat(i, "\n")
    tempParam <- c(coarseOptim[1],refinedOptim[1],planExp3[i,1],refinedOptim[3],planExp3[i,2]) # tempParam <- c(0.02,as.numeric(planExp[1,]))
    temp <- mLogL(optimParam=tempParam,dat=myDb["depthC"],model=myModel,Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
                  NoStat = T,grid=myGrid,db=myDb)
    print(temp)
    planExpLogl3[i] <- temp
  }
  save("planExpLogl3",file='Res/planExpLogl3')
} else{
  load('Res/planExpLogl3')
}

plot(sort(planExpLogl3),type="b")
planExp3[order(planExpLogl3),]
# The possibility to increase r*End values does not improve the likelihood.
# ==> stay with the previous values

myOptimNoStatParamInit = c(coarseOptim[1],refinedOptim)

