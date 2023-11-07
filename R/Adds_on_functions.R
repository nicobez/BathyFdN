source("R/Adds_on_functions/Variance.txt")
source("R/Adds_on_functions/Replication_unfolded_data.txt")
source("R/Adds_on_functions/Projection_point_on_line.txt")
# Modif du 30 05 2022
# unfolding() does everything in a single run with d.across < 0 generated for points that fall in the island
# source("R/Adds_on_functions/Trapeze_membership_number.txt")
source("R/Adds_on_functions/Unfolding.txt")
source("R/Adds_on_functions/Inside.txt")

#############################################
### Vectorial and scalar product in R2
#############################################
f.pv <- function(u,v){
  u[1]*v[2] - u[2]*v[1]
}
f.ps <- function(u,v){
  u[1]*v[1] + u[2]*v[2]
}


#############################################
### Changing coordinates back and forth from (Long, lat) to (Across,Along)
#############################################
f.geo2unfold <- function(db){
  db$setLocators(c('distAcross','distAlong'),ELoc_X())
}

f.unfold2geo <- function(db){
  db$setLocators(c('longitude','latitude'),ELoc_X())
}

#############################################
### Assign anisotropical parameters to the input DB
#############################################
# The range in the ref anisotropical direction (R.across) fluctuates linearly 
#   from r.0 for Dist.across=0 to r.1 = for Dist.across=1
#   R.across(d<=0) = r.0 ; R.across(d>=1) = r.1 = r.0+pente
#    
# The range in the orthogonal direction R.along is R.across/coef.
# The anisotropical coef is assumed constant.

f.aniso <- function(db,r0,r1,coef0,coef1){
  d <- db$getcolumn('distAcross')
  slopeR <- r1-r0
  slopeCoef <- coef1-coef0
  tempR <- (r0+slopeR*d)*(d>=0 & d<=1) + r1*(d>1) + r0*(d<0) # 0*(d<0) # range = 0 for points inside the island ==> no impact of insiders 
  tempCoef <- (coef0+slopeCoef*d)*(d>=0 & d<=1) + coef0*(d<0) + coef1*(d>1)
  db["rAcross"] <- tempR #db$addColumn(tempR,'rAcross')
  db["rAlong"] <- tempR/tempCoef #db$addColumn(tempR/tempCoef,'rAlong')
  db$setLocators(c("dirAcross","rAcross","rAlong"),ELoc_NOSTAT())
}

#################
# Kriging plots
f.image <- function(grid=myGrid,varName,col=rev(myPalette),polyName=NULL,...){
  nx <- grid$getNX(0)
  image(unique(grid['x1']), unique(grid['x2']),
        if(is.null(polyName)) matrix(grid[varName],nx) 
        else matrix(ifelse(grid['Polygon']==1,grid[varName],NA),nx),
        col=col,
        xlab="Longitude",ylab="Latitude",las=1,...)#,
        #xlim=range(grid['x1'][grid['Polygon']==1]),
        #ylim=range(grid['x2'][grid['Polygon']==1]),...)
}


# Function to set changing model's parameters (the sill being fixed externally)
# This allows to not optimizing the sill

setModel = function(optimParam,model,fixParam)
{
  # optimParam = c(nugget,r1,r2)
  # fixParam = c(sill)
  cova = model$getCova(0)
  #cova$setRanges(abs(param[2:3])) 
  # si on entre avec de paramètres < 0 issus d'optim, les modifier modifie la place dans les gradients
  # autant les concserver et mettre des contraintes de positivité dans optim 
  cova$setRanges(optimParam[2:3])
  #cova$setSill(0,0,abs(param[4]))
  cova$setSill(0,0,fixParam[1])
}

#Function to solve a system Q*x = b from the cholesky decomposition of Q (cholQ)

solveFromChol = function(cholQ,b)
{
  b = forwardsolve(t(cholQ),b)
  b = backsolve(cholQ,b)
  b
}

#Function to compute the product of the inverse of Sigma by a vector

invSigmaX = function(invD,Aproj,cholQmat,x) #function(nugget,Aproj,cholQmat,x)
{
  #1/nugget * x - 1/nugget^2 * Aproj %*% solveFromChol(cholQmat,t(Aproj)%*%x)
  invD %*% x - invD %*% Aproj %*% solveFromChol(cholQmat,t(Aproj)%*%invD%*%x) 
  
}


#######################
# Function to compute anisotropical ranges as a function of the across distance
# rStart at the coast line d=0
# rEnd at the sheld edge d=1
setAnisoRange <- function(rStart,rEnd){
  d <- myGrid['distAcross']
  slope <- rEnd-rStart
  res <- (d < 0) * 0 +  
    (d >= 0 & d < 1) * (rStart + slope*d) +
    (d > 1) * rEnd
  res
}

#######################
# Functions to compute minus the loglikelihood in the stationary or the non stationary case
# mLogLC : non stationary sill but fix proportion of nugget ; Q is computed as Dsill %*% Q %*% Dsill
# mLogLSt : stationary sill = 1 but non stationary nugget 

mLogLDepthC=function(optimParam,dat,model,Aproj,mesh,mean = NA,fixParam=1,NoStat=FALSE,grid=NULL,db=NULL)
{
  # grid must be specified and must include NoStat locators in the correct order if NoStat = TRUE
  # This is the user responsability to insure that the grid and the mesh are consistent
  # The nugget is a fix proportion of the non stationary sill
  ndat <- length(dat)
  propNugget = optimParam[1]
  if(!NoStat) setModel(optimParam,model,fixParam)
  if(NoStat){
    r1Start <- optimParam[2]
    r1End <- optimParam[3]
    r2Start <- optimParam[4]
    r2End <- optimParam[5]
    grid['r1'] <- setAnisoRange(r1Start,r1End) 
    grid['r2'] <- setAnisoRange(r2Start,r2End) 
    grid['sill'] <- driftSd(grid['distAcross'])**2
    grid$setLocator("dirAcross",ELoc_NOSTAT(),0)
    grid$setLocator("r1",ELoc_NOSTAT(),1)
    grid$setLocator("r2",ELoc_NOSTAT(),2)
    #grid$setLocator("sill",ELoc_NOSTAT(),3)
    #nostat = NoStatArray(c("A1","R1","R2","V1"),grid)
    nostat = NoStatArray(c("A1","R1","R2"),grid)
    Model_addNoStat(model,nostat)
  }
  
  # building D and invD
  db['sill'] <- driftSd(db['distAcross'])**2
  db['nugget'] <- propNugget * db['sill']
  D <- Diagonal(length(db['nugget']),db['nugget'])
  invD <- Diagonal(length(db['nugget']),1/db['nugget'])
  
  #Precision matrix and Qmat = Q + t(A) %*% invD %*% A
  precOp = PrecisionOpCs(mesh,model)
  Q = cs_toTL(precOp$getQ()) #Precision matrix
  sel <- grid['delatedPolygon']
  selRank <- (1:length(sel))[sel==1]
  temp <- rep(NA,length(sel))
  Dsill <- Diagonal(length(selRank),sqrt(grid['sill'][selRank]))
  Q <- Dsill%*%Q%*%Dsill # precision matrix accounting for a non stationary sill
  
  Qmat = Q+t(Aproj)%*%invD%*%Aproj
  
  #Cholesky decomposition
  cholQmat = chol(Qmat)
  cholQ = chol(Q)
  
  
  #Mean computation (if not provided)
  if(is.na(mean))
  {
    ones = rep(1,ndat)
    invSigmaOnes = invSigmaX(invD,Aproj,cholQmat,ones)
    invSigmaDat = invSigmaX(invD,Aproj,cholQmat,dat)
    mean = sum(ones*invSigmaDat)/sum(ones*invSigmaOnes) # kriging estimation of the stationary mean
    datc = dat -mean
    invSigmaDatc = invSigmaDat-mean*invSigmaOnes
  }else
  {
    datc = dat-mean    
    invSigmaDatc = invSigmaX(invD,Aproj,cholQmat,datc) 
  }
  
  #log det Sigma
  logDetSigma = 2*sum(log(diag(cholQmat)))-2*sum(log(diag(cholQ))) + 2*sum(log(sqrt(diag(D))))  #ndat * log(nugget)
  
  #Quadratic term
  quad = sum(datc*invSigmaDatc)
  
  # minus log likelihood
  0.5 * (quad + logDetSigma)
}



mLogLDepthStd=function(optimParam,dat,model,Aproj,mesh,mean = NA,fixParam=1,NoStat=FALSE,grid=NULL,db=NULL)
{
  # grid must be specified and must include NoStat locators in the correct order if NoStat = TRUE
  # This is the user responsability to insure that the grid and the mesh are consistent
  # the nugget is a fix proportion of the non stationary sill
  ndat <- length(dat)
  if(!NoStat) {
    setModel(optimParam,model,fixParam)
    # building D and invD
    nugget = optimParam[1]
    D <- Diagonal(ndat,nugget)
    invD <- Diagonal(ndat,1/nugget)
  }
  if(NoStat){
    nugStart <- optimParam[1]
    nugEnd <- optimParam[2]
    r1Start <- optimParam[3]
    r1End <- optimParam[4]
    r2Start <- optimParam[5]
    r2End <- optimParam[6]
    grid['r1'] <- setAnisoRange(r1Start,r1End) 
    grid['r2'] <- setAnisoRange(r2Start,r2End) 
    #grid['sill'] <- driftSd(grid['distAcross'])**2
    grid$setLocator("dirAcross",ELoc_NOSTAT(),0)
    grid$setLocator("r1",ELoc_NOSTAT(),1)
    grid$setLocator("r2",ELoc_NOSTAT(),2)
    #grid$setLocator("sill",ELoc_NOSTAT(),3)
    #nostat = NoStatArray(c("A1","R1","R2","V1"),grid)
    nostat = NoStatArray(c("A1","R1","R2"),grid)
    Model_addNoStat(model,nostat)
    
    # building D and invD
    distAcross <- db['distAcross']
    db['nugget'] <- (distAcross < 0) * 0 +  
      (distAcross >= 0 & distAcross < 1) * (nugStart + (nugEnd-nugStart)*distAcross) +
      (distAcross > 1) * nugEnd
    D <- Diagonal(length(db['nugget']),db['nugget'])
    invD <- Diagonal(length(db['nugget']),1/db['nugget'])
  }
  
  #Precision matrix Q and Qmat = Q + t(A)%*% invD %*% A
  precOp = PrecisionOpCs(mesh,model)
  Q = cs_toTL(precOp$getQ()) #Precision matrix
  Qmat = Q+t(Aproj)%*%invD%*%Aproj
  
  #Cholesky decomposition
  cholQmat = chol(Qmat)
  cholQ = chol(Q)
  
  #Mean computation (if not provided)
  if(is.na(mean))
  {
    ones = rep(1,ndat)
    invSigmaOnes = invSigmaX(invD,Aproj,cholQmat,ones)
    invSigmaDat = invSigmaX(invD,Aproj,cholQmat,dat)
    mean = sum(ones*invSigmaDat)/sum(ones*invSigmaOnes) # kriging estimation of the stationary mean
    datc = dat -mean
    invSigmaDatc = invSigmaDat-mean*invSigmaOnes
  }else
  {
    datc = dat-mean    
    invSigmaDatc = invSigmaX(invD,Aproj,cholQmat,datc) 
  }
  
  #log det Sigma
  logDetSigma = 2*sum(log(diag(cholQmat)))-2*sum(log(diag(cholQ))) + 2*sum(log(sqrt(diag(D))))  #ndat * log(nugget)
  
  #Quadratic term
  quad = sum(datc*invSigmaDatc)
  
  # minus log likelihood
  0.5 * (quad + logDetSigma)
}

#####################################
#' Wrapper to optim to save iterations.
#'
#' Source : https://rdrr.io/github/jbryer/visualMLE/src/R/optim_save.R
#' 
#' This function wraps the \link{stats::optim} function and saves the parameters
#' and likelihood estimation at each step of the algorithm.
#'
#' @param par initial parameters to to be optimized over.
#' @param fn the function to minimized (or maximized).
#' @param ... other parameters passed to optim.
#' @return the results of optim with two additional elements, iterations with a
#'         a list of the values at each iteration of the algorithm and
#'         iterations_df which is a data.frame version of the list.
#' @seealso stats::optim
#' @export
optim_save <- function(par, fn, ...) {
  iterations <- list()
  wrap_fun <- function(parameters, ...) {
    n <- length(iterations)
    result <- fn(parameters, ...)
    iterations[[n + 1]] <<- c(parameters, result)
    return(result)
  }
  optim_out <- stats::optim(par, wrap_fun, method = c("Nelder-Mead"), ...)
  optim_out$iterations <- iterations
  optim_out$iterations_df <- as.data.frame(do.call(rbind, iterations))
  names(optim_out$iterations_df) <- c(paste0('Param', 1:length(par)), 'Result')
  optim_out$iterations_df$Iteration <- 1:nrow(optim_out$iterations_df)
  return(optim_out)
}

