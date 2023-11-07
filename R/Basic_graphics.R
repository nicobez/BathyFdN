myDb <- depthRegul10000$clone()

# Plot data in geographical space ---------------------------------------------------------------
f.unfold2geo(myDb)
myDb$setLocator('depth',ELoc_Z(),0)

tempSegment <- as.data.frame(list(x=ref.line$x,y=ref.line$y,xend=ref.line.2$x,yend=ref.line.2$y))

tempDf <- as.data.frame(list(x=myDb['longitude'],y=ref.line$y,xend=ref.line.2$x,yend=ref.line.2$y))
tempGraph <- plot.point(myDb,name_color='depth',name_size = "depth",
                        legend.name.color = "depth",
                        show.legend.symbol =T)#,mode=2,xlab='Longitude',ylab='Latitude'
                        
tempGraph <- tempGraph +
  theme_bw()+
  scale_color_gradientn(colours=myPalette)+
  geom_polygon(data = as.data.frame(ref.line),aes(x=x,y=y),fill="grey",color=1)+
  geom_path(data = as.data.frame(ref.line.2),aes(x=x,y=y),color=1,linewidth=1)+
  geom_segment(data=tempSegment,aes(x=x,y=y,yend=yend,xend=xend),linewidth=1)

print(tempGraph)

dev.print(device = png, file = "Res/regul.data.png",width=450,height=450)

# Plot data in the unfolded space ---------------------------------------------------------------
f.geo2unfold(myDb)

myDf <- cbind(myDb['distAcross'],myDb['distAlong'],myDb['depth'])

tempGraph <- plot.point(myDb,name_color='depth',name_size = "depth",
                        legend.name.color = "depth",
                        show.legend.symbol =T)
tempGraph <- tempGraph +
  theme_bw()+
  scale_color_gradientn(colours=myPalette)

print(tempGraph)
dev.print(device = png, file = "Res/regul.data.unfolded.png",width=450,height=450)

rm(tempN,tempGraph,tempSegment,myDf)

myDb <- depthRegul$clone()


