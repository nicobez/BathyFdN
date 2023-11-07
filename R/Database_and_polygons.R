
# # Reading polygons and reference lines ------------------------------------
# NEEDS RGeostats !!!!!!!!!!!!!!!!!
# # lecture du polygone 
# myPolygonShelf <- readRDS("Data/my.polygone_shelf.rds")@sets[[1]]
# 
# # creating de 2 polygons : 1 for FNA ;  1 for the shelfbreak.
# temp <- Db_create()
# temp$addColumns(tab=myPolygonShelf$x[c(1:25,48,49)],"Long",ELoc_X(),0)
# temp$addColumns(tab=myPolygonShelf$y[c(1:25,48,49)],"Lat",ELoc_X(),1)
# polyFernando <- Polygons_createFromDb(temp)
# 
# temp <- Db_create()
# temp$addColumns(tab=myPolygonShelf$x[26:47],"Long",ELoc_X(),0)
# temp$addColumns(tab=myPolygonShelf$y[26:47],"Lat",ELoc_X(),1)
# polyShelf <- Polygons_createFromDb(temp)

source("Data/ref.line.R")
source("Data/ref.line.2.R")
load("Data/kriging.polygon.RData")

temp <- Db_create()
temp$addColumns(tab=res$x,"Long",ELoc_X(),0)
temp$addColumns(tab=res$y,"Lat",ELoc_X(),1)
polyKriging <- Polygons_createFromDb(temp)
  
rm(temp)

# Reading 1 ping data ------------------------------------------------------------
print("Reading data")
df <- read.csv(file = 'Data/FAROFA123_logSa_1ping.csv')
df = na.omit(df)

# Cleaning dataset : removing 5 isolated points
sel <- df$Longitude > -32.41892 & df$Longitude < -32.41783 & df$Latitude > -3.799620 & df$Latitude < -3.798449
df <- df[!sel,]
names(df) <- c("longitude","latitude","sa70","sa200","depth","time","esdu","farofa")

# creation DB
depthHighDef <- Db_create()
depthHighDef$addColumns(df$longitude,'longitude',ELoc_X(),0)
depthHighDef$addColumns(df$latitude,'latitude',ELoc_X(),1)
depthHighDef$addColumns(df$depth,'depth',ELoc_Z(),0)

# Regularization ----------------------------------------------------------
# print("Regularization")
# pixel.size <- 0.0001 # after regularization it remains 82523 data points
# xmin <- c(min(db[,db.getcols(db,"x",1)]),min(db[,db.getcols(db,"x",2)]))
# xmax <- c(max(db[,db.getcols(db,"x",1)]),max(db[,db.getcols(db,"x",2)]))
# dx1 <- xmax[1]-xmin[1]
# dx2 <- xmax[2]-xmin[2]
# marge <- 0.01
# xmin = xmin - marge * c(dx1,dx2)
# xmax = xmax + marge * c(dx1,dx2)
# dx1 <- xmax[1]-xmin[1]
# dx2 <- xmax[2]-xmin[2]
# regul.grid <- db.create(flag.grid=T,
#                         x0=xmin,
#                         dx=c(pixel.size,pixel.size),
#                         nx=c(round(dx1/pixel.size,0),round(dx2/pixel.size,0)))
# #regul.grid <- db.grid.init(db,dcell=pixel.size)
# regul.grid <- db.stat.grid(db,regul.grid,names=4:db$natt,#names=c("Depth","FAROFA"),
#                            radix="",modify.target = F)
# regul.grid <- db.rename(regul.grid,c("x1","x2"),c("Longitude","Latitude"))
# 
# # Formatting regularized data as set of points (and not as gridded points)
# db.regul <- db.grid2point.sampling(regul.grid,npack=1,names=4:regul.grid$natt)#c("Depth","FAROFA"))
# db.regul <- db.reduce(db.sel(db.regul,!is.na(Depth)))
# db.regul <- db.delete(db.regul,'sel')
# db.regul <- db.locate(db.regul,"Depth","z")
# 
# df.regul <- db.regul@items[,2:db.regul$natt]
# 
# save(df.regul,file="Data/df.regul.RData")

print("Loading regularized data")
load("Data/df.regul.RData")

names(df.regul) <- c("longitude","latitude","sa70","sa200","depth","time","esdu","farofa")

# creation DB
depthRegul <- Db_create()
depthRegul$addColumns(df.regul$longitude,'longitude',ELoc_X(),0)
depthRegul$addColumns(df.regul$latitude,'latitude',ELoc_X(),1)
depthRegul$addColumns(df.regul$depth,'depth',ELoc_Z(),0)

# Unfolding 1 ping data --------------------------------------------------------------
# By convention Dist.across ranges from 0 to 1
# If required it could be rescaled in order to get roughly similar scales along et across
# res.proj$Dist.across <- res.proj$Dist.across/20 

# ### If required, projections can be recomputed (~30 mn)
# res.proj <- f.unfolding(db$Longitude,db$Latitude,ref.line = ref.line,ref.line.2=ref.line.2)
# res.proj <- f.unfolding(db.10000$Longitude,db.10000$Latitude,ref.line = ref.line,ref.line.2=ref.line.2)
# save(res.proj,file="Data/proj.1ping.RData")

### Otherwise they can be loaded ==> res.proj and res.trapeze 
load("Data/proj.1ping.RData") 

depthHighDef$addColumns(res.proj$Dist.across,'distAcross')
depthHighDef$addColumns(res.proj$Dist.along,'distAlong')
depthHighDef$addColumns(res.proj$Dir.across,'dirAcross')
depthHighDef$addColumns(res.proj$Trapeze.nb,'trapezeNb')

# Unfolding regularized data --------------------------------------------------------------
### If required, projections can be recomputed (~ 2 mn)
# res.proj.regul <- f.unfolding(db.regul$Longitude,db.regul$Latitude,ref.line = ref.line,ref.line.2=ref.line.2)
# save(res.proj.regul,file="Data/proj.regul.data.RData")

### Otherwise they can be loaded  
load("Data/proj.regul.data.RData") 

depthRegul$addColumns(res.proj.regul$Dist.across,'distAcross')
depthRegul$addColumns(res.proj.regul$Dist.along,'distAlong')
depthRegul$addColumns(res.proj.regul$Dir.across,'dirAcross')
depthRegul$addColumns(res.proj.regul$Trapeze.nb,'trapezeNb')

# # Esdu 25 m ----------------------------------------------------------
# # Be careful it does not come with unfolded coordinates
# df <- read.csv('Data/mel_nig_F123_ESDU_25m.csv',header = T)
# db.esdu25 <- db.create(df)
# 
# db.esdu25 <- db.rename(db.esdu25,c("Horizontal_sa_surface_fish_70","Bottom"),c("Sa70","Depth"))
# db.esdu25 <- db.locate(db.esdu25,c(2,3),loctype="x")
# db.esdu25 <- db.locate(db.esdu25,c("Depth"),"z")

########### 
myPolygon = Polygons()
polyset1 = PolySet(x=ref.line$x, y = ref.line$y)
polyset2 = PolySet(x=ref.line.2$x, y = ref.line.2$y)

myPolygon$addPolySet(polyset1)
myPolygon$addPolySet(polyset2)


