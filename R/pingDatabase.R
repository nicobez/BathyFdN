######################################
# Read the ping Dataset after they have been cleaned.
# The coordinates lack precision ==> lost of duplicates.
# ==> regularization towards small pixels : build a grid + mean per cell + turn the output as a point Db by removing all the NAs.
# Then, unfold the regularized dataset.
#####################################
df <- read.table("Data/FAROFA123_logSa_1ping_cleaned",header=T)
# creation DB
depthPing <- Db_create()
depthPing$addColumns(df$longitude,'longitude',ELoc_X(),0)
depthPing$addColumns(df$latitude,'latitude',ELoc_X(),1)
depthPing$addColumns(df$depth,'depth',ELoc_Z(),0)
depthPing$addColumns(df$farofa,'farofa')

# GRID 
pixelSize <- 0.0001 
marge <- 0.1 # % of extension
xmin <- c(min(df$longitude),min(df$latitude)) #c(-32.55,-3.93) #
xmax <- c(max(df$longitude),max(df$latitude)) #c(-32.3,-3.77) #
dx1 <- xmax[1]-xmin[1]
dx2 <- xmax[2]-xmin[2]
xmin = xmin - marge * c(dx1,dx2)
xmax = xmax + marge * c(dx1,dx2)
dx1 <- xmax[1]-xmin[1]
dx2 <- xmax[2]-xmin[2]
nx <- round(dx1/pixelSize,0)
ny <- round(dx2/pixelSize,0)
nxy <- c(nx,ny)
myGrid <- DbGrid_create(x0=xmin,dx=c(pixelSize,pixelSize),nx=nxy)

# Regularisation depth
temp <- dbStatisticsPerCell(depthPing,myGrid,oper=EStatOption(1),name1='depth',name2='depth')
tempSel <- !is.na(temp)
depthRegul <- Db_create()
depthRegul$addColumns(myGrid['x1'][tempSel],'longitude',ELoc_X(),0)
depthRegul$addColumns(myGrid['x2'][tempSel],'latitude',ELoc_X(),1)
depthRegul$addColumns(temp[tempSel],'depth',ELoc_Z(),0)

# Regularisation farofa
temp <- dbStatisticsPerCell(depthPing,myGrid,oper=EStatOption(1),name1='farofa',name2='farofa')
tempSel <- !is.na(temp)
depthRegul$addColumns(round(temp[tempSel],0),'farofa')

ggplot() + plot(depthRegul,name_color="depth")

# Unfolding regularized data (~2 mn)--------------------------------------------------------------
# By convention Dist.across ranges from 0 to 1
# If required it could be rescaled in order to get roughly similar scales along et across
# res.proj$Dist.across <- res.proj$Dist.across/20 
start <- Sys.time()
depthRegulProj <- f.unfolding(depthRegul['longitude'],depthRegul['latitude'],ref.line = ref.line,ref.line.2=ref.line.2)
write.table(depthRegulProj,"Data/depthRegulProj",row.names = F)
print( Sys.time() - start )
#depthRegulProj <- read.table("Data/depthRegulProj",header=T)

depthRegul$addColumns(depthRegulProj$Dist.across,'distAcross')
depthRegul$addColumns(depthRegulProj$Dist.along,'distAlong')
depthRegul$addColumns(depthRegulProj$Dir.across,'dirAcross')
depthRegul$addColumns(depthRegulProj$Trapeze.nb,'trapezeNb')

# Saving output
varNames <- c('longitude','latitude','depth','farofa','distAcross','distAlong','dirAcross','trapezeNb')
temp <- as.data.frame(depthRegul$getColumnsAsVVD(varNames),col.names = varNames)
write.table(temp,"Data/FAROFA123_depth_cleaned_regul",row.names=F)


