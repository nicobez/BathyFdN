
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
polyTotal <- read.table('Data/polygonComplet.ascii',header=T)

########### Create a polygon based on the lines contouring all the data

myPoly = Polygons()
polyset1 = PolySet(x=polyTotal$x,y=polyTotal$y)
myPoly$addPolySet(polyset1)

# 
# temp <- Db_create()
# temp$addColumns(tab=polyTotal$x,"Long",ELoc_X(),0)
# temp$addColumns(tab=polyTotal$y,"Lat",ELoc_X(),1)
# myPoly <- Polygons_createFromDb(temp) ### CAREFUL ### : this takes the convex hull of the data/polygon
# rm(temp)

########### Create a polygon based on the reference lines used in the unfolding operation
myPolyProj = Polygons()
polyset1 = PolySet(x=ref.line$x, y = ref.line$y)
polyset2 = PolySet(x=ref.line.2$x, y = ref.line.2$y)
myPolyProj$addPolySet(polyset1)
myPolyProj$addPolySet(polyset2)

