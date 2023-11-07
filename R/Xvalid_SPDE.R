### Loading polygon
xvalid.polygon <- polygon.read('Data/xvalid.polygon.ascii')

z.in <- db.extract(db.polygon(wdb,polygon=xvalid.polygon),names='Depth.residuals',flag.compress = T)

coef.multiplicatif <- c(0.25,0.5,seq(1,3,0.5))

res.xvalid.transect.spde <- array(NA,dim=c(length(z.in),length(coef.multiplicatif)))

for(k in 1:length(coef.multiplicatif)){ # i <- 1.5
    i <- coef.multiplicatif[k]
    cat("i = ",i,"\n")
    wdb.2 <- db.add(wdb,Range1=i*(r1.deb+Dist.across*20*(r1.fin-r1.deb)))
    wdb.2 <- db.add(wdb.2,Range2=i*(r2.deb+Dist.across*20*(r2.fin-r2.deb)))
    wdb.2 <- db.add(wdb.2,Dir.aniso=Dir.across)
    wdb.2 <- db.locate(wdb.2,c("Dir.aniso","Range1","Range2"),"nostat")
    
    spde.grid <- db.add(grid,Range1=i*(r1.deb+Dist.across*20*(r1.fin-r1.deb)))
    spde.grid <- db.add(spde.grid,Range2=i*(r2.deb+Dist.across*20*(r2.fin-r2.deb)))
    spde.grid <- db.add(spde.grid,Dir.aniso=Dir.across)
    spde.grid <- db.locate(spde.grid,c("Dir.aniso","Range1","Range2"),"nostat")
    
    spde.grid.vect <- db.create(db.extract(spde.grid,names=c('x1','x2',"Range1","Range2","Dir.aniso")))
    spde.grid.vect <- db.locate(spde.grid.vect,c("Dir.aniso","Range1","Range2"),"nostat")
    
    wdb.in <- db.reduce(db.polygon(wdb.2,polygon=xvalid.polygon))
    wdb.out <- db.reduce(db.polygon(wdb.2,polygon=xvalid.polygon,flag.out=T))
    
    ### SPDE xvalid
    wdb.in <- db.locate(wdb.in,c("Dir.aniso","Range1","Range2"),"nostat")
    tmp <- db.extract(wdb.in,names=c('Longitude','Latitude',"Dir.aniso","Range1","Range2"))
    names(tmp) <- c('x1','x2',"Dir.aniso","Range1","Range2")
    tmp <- db.append(spde.grid.vect,tmp)
    
    Q= spde.matrices(,tmp,model.vario.spde,flag.Q = T,triswitch = "Q",nostat=c("M1A","M1R1","M1R2"))$Q 
    me=meshing(,tmp,triswitch = "Q")#,gext=0.1)
    B=mesh.barycenter(wdb.out,me,flag.sparse = T)
    sol=Matrix:::solve(nug*Q+t(B)%*%B,t(B)%*%wdb.out$Depth.residuals)[,1]
    
    res.xvalid.transect.spde[,k] <- sol[spde.grid.vect$nech+1:wdb.in$nech]
}

plot(z.in,res.xvalid.transect.spde[,1],pch=20,cex=2,cex.axis=2,
     xlim=range(c(z.in,res.xvalid.transect.spde)),ylim=range(c(z.in,res.xvalid.transect.spde)),
     #xlab="True value",ylab="Re-estimated value",
     xlab="",ylab="",las=1)
for(i in 2:length(coef.multiplicatif)) points(z.in,res.xvalid.transect.spde[,i],pch=20,cex=2,col=i)

abline(0,1,lty=1,lwd=1)
abline(v=mean(z.in),h=mean(z.in),lty=1,lwd=2,col=1)

dev.print(device = png, file = "Res/xvalid.png",width=800,height=800)


