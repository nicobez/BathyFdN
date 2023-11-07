#### GRID
pixelSize <- 0.001 #0.0005 #.001 
marge <- 0.15 # % of extension

# xmin <- c(min(myDb['longitude']),min(myDb['latitude'])) #c(-32.55,-3.93) #
# xmax <- c(max(myDb['longitude']),max(myDb['latitude'])) #c(-32.3,-3.77) #
xmin <- c(min(polyTotal$x),min(polyTotal$y))
xmax <- c(max(polyTotal$x),max(polyTotal$y))
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

# Adding unfolded coordinates (~ 40 secondes )
temp <- f.unfolding(myGrid['x1'],myGrid['x2'],ref.line,ref.line.2)

myGrid['distAlong'] <- temp$Dist.along
myGrid['distAcross'] <- temp$Dist.across
myGrid['dirAcross'] <- temp$Dir.across

# The polygon is delated to handle borders effects of SPDE
# myGridExtPolygon <- myGridExt$clone()

db_polygon(myGrid, myPoly)

sel <- myGrid['Polygon']
# IN = GRAIN = 1 #11258 ; OUT = PORE = 0 #37142
image <- BImage(nxy)
imageRes <- BImage(nxy)
tempVD <- VectorDouble(nx*ny)
morpho_double2image(nxy,sel,vmin=1,vmax=2,imagout=image)
morpho_dilation(1, c(20,20), image, imageRes)
morpho_image2double(imageRes, 0, 1., 0., tempVD)
myGrid$addColumns(tempVD,"delatedPolygon",ELoc_SEL(),0)

# Figure
f.image(varName='delatedPolygon',col=c("white",myPalette[7]))
polygon(polyTotal,lwd=2,col="lightblue")
points(myDb['longitude'][myDb['farofa']<3],myDb['latitude'][myDb['farofa']<3],pch=20,cex=0.05,col='lightcoral')
points(myDb['longitude'][myDb['farofa']==3],myDb['latitude'][myDb['farofa']==3],pch=20,cex=0.05,col='lightblue4')
lines(ref.line.2,col=1,lwd=2)
polygon(ref.line,col="white",border = 1,lwd=2)
for(i in 1:length(ref.line.2$x)) segments(ref.line$x[i],ref.line$y[i],ref.line.2$x[i],ref.line.2$y[i],lty = 1,col=1,lwd=2)

for(j in seq(0.1,1.2,0.1)){
  for(i in 1:length(ref.line.2$x)) {
 x0 <- ref.line$x[i] + (ref.line.2$x[i] - ref.line$x[i])*j
 y0 <- ref.line$y[i] + (ref.line.2$y[i] - ref.line$y[i])*j
 x1 <- ref.line$x[i+1] + (ref.line.2$x[i+1] - ref.line$x[i+1])*j
 y1 <- ref.line$y[i+1] + (ref.line.2$y[i+1] - ref.line$y[i+1])*j
 segments(x0,y0,x1,y1,lty = 2,col=1,lwd=1)
  }
}

# Important remark : 
# when computing omni directional variograms by slices of distAcross,
# g(h) for h > 0.005 = 5*pixelSize  ==> along distance structure up to "corde effect"

# Meshing
myMesh <- MeshETurbo(myGrid)




