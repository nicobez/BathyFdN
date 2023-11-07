# Stationary SPDE 
# First with user's defined parameters and then with nugget and ranges parameters estimated by ML (the sill is fixed = 1)
# This is justified as we are working on the standardized depth residuals

#######################
# Stationary KO with user's defined parameters
nugget = 0.05
r1 = 0.01
r2 = 0.02
sill = 1.
myParam = c(nugget,r1,r2,sill)

myModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=myParam[c(2,3)],sill=myParam[4])

# Direct calculus 
Q = cs_toTL(PrecisionOpCs(myMesh,myModel)$getQ()) #Precision matrix
M1 <- nugget^2*Q + t(Aproj)%*%Aproj
M2 <- t(Aproj)%*%myDb["depthResiduals"]
zKri <-  Matrix:::solve(M1,M2)
myGridExt['kriStat'] <- zKri*driftSd(myGridExt['distAcross'])+driftM(myGridExt['distAcross'])

p = ggplot()
p = p + plot.grid(myGridExt, name_raster='kriStat',palette='viridis',limits=c(0,100))
ggPrint(p)

#######################
# Stationary KO with ML parameters
# The sill is fixed to 1 and is NOT optimized

myOptimParamInit = c(nugget,r1,r2)
start <- Sys.time()
mLogL0 = abs(mLogL(myOptimParamInit,myDb["depthResiduals"],myModel,Aproj,myMesh))
print( Sys.time() - start )
# 6.8 seconds to compute a logL

# Defining the cost function to be optimzed by optim()
# The logL is normalised by the initial loglikelihood (to prevent numerical issues)
fOpt = function(optimParam,dat,model,Aproj,mesh,mean=NA,fixParam)
{
  print(optimParam)
  res = mLogL(optimParam,dat,model,Aproj,mesh,mean,fixParam)/mLogL0 
  print(res)
  return(res)
}

myMaxit <- 10
start <- Sys.time()
resOptimStat = optim_save(myOptimParamInit,fOpt,model=myModel,dat = myDb["depthResiduals"],Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
            lower=rep(0.001,4),upper=rep(0.5,4),control = list(maxit=myMaxit))
print( Sys.time() - start )

essai <- optim(myOptimParamInit,mLogL,model=myModel,dat = myDb["depthResiduals"],Aproj=Aproj,mesh=myMesh,mean=0,fixParam=1,
           lower=rep(0.001,4),upper=rep(0.5,4),control = list(maxit=1))



# The number of access to fOpt is much larger than maxit.
# The result corresponds to the output dataframe at the beginning of each batch of 7 accesses to fOpt.
# The current value is the one at the beginning of each batch.
# The following 6 accesses to fOpt are probably there to prepare the next step (gradient, ...).
# The optimized parameters are thus those available at the beginning of the final batch of computations.
# 39 minutes for maxit = 10 , i.e. around the elementary computing time tiumes 400 
# 2 hours of computation without maxit limitation
#
# Initial values
# 0.05 0.01 0.02
# Optimized values
# 0.02113952 0.01666842 0.01816150
# r1 ~ r2 ==> isotropy

# save("resOptimStat",file='Res/resOptimStat')
load('Res/resOptimStat')

plot(resOptimStat$iterations_df$Param1[c(1,50,57,71+7*(0:(myMaxit-4)))],type="l")
par(new=T)
plot(resOptimStat$iterations_df$Param2[c(1,50,57,71,78)],col=2,type="l")

lines(resOptimStat$iterations_df$Param3[c(1,50,57,71,84)],col=3)
plot(resOptimStat$iterations_df$Result[c(1,50,57,71,84)],type="l")

myParamStat <- c(mean(resOptimStat$iterations_df$Param1[100:196]),
                 mean(resOptimStat$iterations_df$Param2[100:196]),
                 mean(resOptimStat$iterations_df$Param3[100:196]))
myParamStat <- resOptimStat$values


######################

myModel = Model_createFromParam(type=ECov_BESSEL_K(),param=1,ranges=myParamStat[c(2,3)],sill=1)

Q = cs_toTL(PrecisionOpCs(myMesh,myModel)$getQ()) #Precision matrix
M1 <- myParamStat[1]^2*Q + t(Aproj)%*%Aproj
M2 <- t(Aproj)%*%myDb["depthResiduals"]
zKri <-  Matrix:::solve(M1,M2)
myGridExt['kriStatMaxLike'] <- zKri*driftSd(myGridExt['distAcross'])+driftM(myGridExt['distAcross'])

p = ggplot() + plot.grid(myGridExt, name_raster='kriStatMaxLike',palette='viridis',limits=c(0,100))
ggPrint(p)
