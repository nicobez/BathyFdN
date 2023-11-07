
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
# Coarse experimental design 
#############################
if(compute){
  nugget <- c(0.01,0.05,0.1)
  r1 <- c(0.01,0.025,0.05)
  r2 <- c(0.01,0.025,0.05)
  planExpStatStd1 <- expand.grid(nugget,r1,r2)
  planExpStatStd1[,4] <- rep(NA,NROW(planExpStatStd1))
  names(planExpStatStd1) <- c("nugget","r1","r2","mLogL")
  for(i in 1:NROW(planExpStatStd1)){
    cat("\n", i, "/", NROW(planExpStatStd1), "\n")
    temp <- mLogLDepthStd(optimParam=as.numeric(planExpStatStd1[i,1:3]),dat=myDb["depthStd"],
                          model=myModel,Aproj=Aproj,mesh=myMesh,
                          mean=0,fixParam=1,NoStat = F,grid=myGrid,db=myDb)
    planExpStatStd1[i,4] <- temp
  }
  planExpStatStd1 <- planExpStatStd1[order(planExpStatStd1[,4]),]
  save("planExpStatStd1",file='Res/planExpStatStd1')
} else{
  load('Res/planExpStatStd1')
}

plot(planExpStatStd1[,4],type="b")

# Comments:
# the minimum LogL correspond all to nug = 0.01 then 0.05 and then 0.1
# Then r1 = r2 = 0.025 are the prefered range values for each possible nugget.
# The output is robust anb OK (no need to refine the design)

myOptimStatParamStd = as.numeric(planExpStatStd1[1,1:3])

