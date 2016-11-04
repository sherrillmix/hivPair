source('prinComp.R')

selector<-!apply(is.na(hiv[,names(goodTargetCols)]),1,any)
tmp<-hiv[selector,names(goodTargetCols)]
rownames(tmp)<-hiv[selector,'Renamed']
pca<-prcomp(tmp,scale.=TRUE)
select<-1:3
pcaPoints<-t(t(pca$x)/pca$sdev/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
pointLim <- range(-pcaPoints[,select])
pcaArrows<-t(t(-pca$rotation[,select])*pca$sdev[select]*sqrt(nrow(pca$x)))	#figure out the arrow positions based on loadings scaled by sdev
arrowLim<-range(pcaArrows)
ratio <- max(arrowLim/xlim)

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
apply(pcaArrows,1,function(xx)arrow3d(c(0,0,0),xx/ratio*.8,type='rotation',col='red',n=50,width=.1,thickness=.1,barblen=.02))
text3d(pcaArrows[,1]/ratio*.8,pcaArrows[,2]/ratio*.8,pcaArrows[,3]/ratio*.8,sub('\\(.*$','',targetCols[rownames(pcaArrows)]))
for(ii in 1:90){
  message(ii)
  view3d(userMatrix=rotationMatrix(2*pi*ii/90,1,-1,-1))
  rgl.snapshot(filename=sprintf('out/pcaAnim/%03d.png',ii))
}
system('convert -delay 8 -loop 0 -layers optimize out/pcaAnim/*.png out/pcaAnim/pca.gif')
