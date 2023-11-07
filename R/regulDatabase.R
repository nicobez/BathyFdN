
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
depthRegul$addColumns(df.regul$farofa,'farofa')

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

# Cleaning truncated data due to echosounder's threshold 
tempZ <- unique(df.regul$depth)
tempFreq <- sapply(1:length(tempZ), function(i) sum(df.regul$depth==tempZ[i]))
summary(tempFreq) # More than 50% of the observed values are observed only once
plot(tempZ,tempFreq,ylim=c(0,50))
tempZSel <-  (tempZ > 190 & tempZ < 210 & tempFreq > 50) | 
  (tempZ > 290 & tempZ < 310 & tempFreq > 50) |
  tempZ > 390

tempZKept <- tempZ[!tempZSel]
tempFreqKept <- tempFreq[!tempZSel]
plot(tempZKept,tempFreqKept)

tempZToBeRemoved <- tempZ[tempZSel]

sel <- df.regul$depth %in% tempZToBeRemoved
df.regul <- df.regul[!sel,]
write.table(df.regul,"Data/proj.regul.data.cleaned",row.names=F)
res.proj.regul <- read.table("Data/proj.regul.data.cleaned",header=T)


