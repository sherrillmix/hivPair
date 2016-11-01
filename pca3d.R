source('prinComp.R')


library(scatterplot3d)
cols4<-rainbow.lab(length(unique(hiv$fluidSelectDonor)),alpha=.2)
names(cols4)<-unique(hiv$fluidSelectDonor)
pdf('out/pca3d.pdf')
  s3d<-scatterplot3d(pcaPoints[,1],pcaPoints[,2],pcaPoints[,3],type='h',pch=21,bg=cols[as.character(hiv[selector,'fluidSelectDonor'])],color=cols4[as.character(hiv[selector,'fluidSelectDonor'])],angle=50)
  #zs<-c(0,.2)
  #sapply(zs,function(z)s3d$plane3d(z,0,0))
dev.off()

library('rgl')
cols4<-rainbow.lab(length(unique(hiv$fluidSelectDonor)))
names(cols4)<-unique(hiv$fluidSelectDonor)
plot3d(pcaPoints[,1:3],col=cols4[as.character(hiv[selector,'fluidSelectDonor'])],type='s',radius=.005,box=FALSE,axes=FALSE,xlab='',ylab='',zlab='')
for(ii in 1:180){
  message(ii)
  view3d(userMatrix=rotationMatrix(2*pi*ii/180,1,-1,-1))
  rgl.snapshot(filename=sprintf('out/pcaAnim/%03d.png',ii))
}

