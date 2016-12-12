library(vipor)
library(cluster)
if(!exists('hiv'))source('readData.R')
expandedHull<-function(xys,magnification=1,type=c('convex','ellipse')){
  type<-match.arg(type)
  centroid<-apply(xys,2,mean)
  magnified<-t(apply(xys,1,function(xy)(xy-centroid)*magnification+centroid))
  if(type=='convex')return(magnified[chull(magnified),])
  else if(type=='ellipse')return(predict(ellipsoidhull(magnified)))
  else stop("Unknown type")
}

selector<-!apply(is.na(hiv[,goodTargetCols]),1,any)
xFlip<--1
yFlip<- -1
subselect<-hiv[selector,]$Selection=='UT'
scaled<-apply(hiv[selector,goodTargetCols],2,function(x,subselect)(x-mean(x[subselect]))/sd(x[subselect]),rep(TRUE,sum(selector)))
pca<-prcomp(scaled,scale.=FALSE)
pcaPoints<-t(t(scaled %*% pca$rotation)/pca$sdev/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
importance<-summary(pca)$importance[2,]
cols<-rainbow(length(unique(hiv$fluidSelectDonor)),alpha=.6)
cols2<-rainbow(length(unique(hiv$fluidSelectDonor)),alpha=.8)
cols3<-rainbow(length(unique(hiv$fluidSelectDonor)),alpha=.02)
names(cols)<-names(cols2)<-names(cols3)<-sort(unique(hiv$fluidSelectDonor))
pdf('out/pca.pdf',width=5,height=5)
  for(subfig2 in c(FALSE,TRUE)){
  select<-1:2
  xlim <- ylim <- range(yFlip*pcaPoints[,select])
  pcaArrows<-t(t(pca$rotation[,select])*pca$sdev[select]*sqrt(nrow(pca$x)))	#figure out the arrow positions based on loadings scaled by sdev
  xlimArrow<-range(xFlip*pcaArrows[,1])
  ylimArrow<-range(yFlip*pcaArrows[,2])
  ratio <- max(xlimArrow/xlim, ylimArrow/ylim)
  par(mar=c(3.5,3.5,.1,.1))
  plot(
    1,1,type='n',las=1,
    #xlim=range(-pcaPoints[,select[1]]),ylim=range(-pcaPoints[,select[2]]),
    xlim=xlim*1.1,ylim=ylim*1.1,mgp=c(2.4,.7,0),
    xlab=sprintf('Principal component %d (%d%% of variance)',select[1],round(importance[select[1]]*100)),
    ylab=sprintf('Principal component %d (%d%% of variance)',select[2],round(importance[select[2]]*100))
  )
  if(!subfig2){
    arrows(0,0,xFlip*pcaArrows[,1]/ratio,yFlip*pcaArrows[,2]/ratio,length=.1) #draw arrows
    arrowText<-dimnames(pcaArrows)[[1]] #get rownames of loadings for labels
    par(lheight=.7)
    text(xFlip*pcaArrows[,1]/ratio,yFlip*pcaArrows[,2]/ratio,sub(' ','\n',sub('\\(.*$','',rownames(pcaArrows))),cex=.75) #label the arrows
  }
  donorFluidBinary<-list('Donor Genital Untreated'=c(TRUE,TRUE),'Donor Plasma Untreated'=c(TRUE,FALSE),'Recipient Plasma Untreated'=c(FALSE,FALSE))
  for(dfName in names(donorFluidBinary)){
    donorFluid<-donorFluidBinary[[dfName]]
    donorFluidSelector<-hiv[selector,'isDonor']==donorFluid[1]&hiv[selector,'isGenital']==donorFluid[2]&hiv[selector,'Selection']=='UT'
    hull<-expandedHull(pcaPoints[donorFluidSelector,select],1,'ellipse')
    polygon(xFlip*hull[,1],yFlip*hull[,2],border=cols2[dfName],col=cols3[dfName],lwd=2.4)
  }
  if(subfig2)thisSelect<-!subselect
  else thisSelect<-subselect
  points(xFlip*pcaPoints[,select[1]][thisSelect],yFlip*pcaPoints[,select[2]][thisSelect],bg=cols[as.character(hiv[selector,'fluidSelectDonor'])[thisSelect]],col=cols2[as.character(hiv[selector,'fluidSelectDonor'])[thisSelect]],pch=21,cex=1.5)
  if(subfig2) legend('topright',names(cols),pch=21,col=cols2,pt.bg=cols,inset=.01,pt.cex=1.5,cex=.7)
  }
dev.off()



#distance calculating
scaled<-pcaPoints[,1:2]
sampleNames<-paste(hiv[selector,'sampleFluidSelect'],ave(hiv[selector,'sampleFluidSelect'],hiv[selector,'sampleFluidSelect'],FUN=function(x)1:length(x)))
rownames(scaled)<-sampleNames
dists<-dist(scaled)
distMat<-as.matrix(dists)


xx<-hiv[selector,]
ids<-list(
  donor=sampleNames[xx$isDonor&xx$Selection=='UT'&xx$Fluid=='PL']
  ,donorGenital=sampleNames[xx$isDonor&xx$Selection=='UT'&xx$Fluid!='PL']
  ,donorAlpha=sampleNames[xx$isDonor&xx$Selection=='A2'&xx$Fluid=='PL']
  ,donorBeta=sampleNames[xx$isDonor&xx$Selection=='BE'&xx$Fluid=='PL']
  ,recAlpha=sampleNames[!xx$isDonor&xx$Selection=='A2'&xx$Fluid=='PL']
  ,recBeta=sampleNames[!xx$isDonor&xx$Selection=='BE'&xx$Fluid=='PL']
  ,recipient=sampleNames[!xx$isDonor&xx$Selection=='UT']
)

pairRecDist<-unlist(ave(split(cbind(scaled,!hiv[selector,'isDonor']&hiv[selector,'Fluid']=='PL'&hiv[selector,'Selection']=='UT'),1:nrow(scaled)),hiv[selector,'Pair ID'],FUN=function(df){
  df<-do.call(rbind,df)
  select<-df[,ncol(df)]==1
  df<-df[,-ncol(df)]
  recipientCentroid<-apply(df[select,],2,mean)
  recipientDiff<-t(t(df)-recipientCentroid)
  return(apply(recipientDiff^2,1,sum))
}))

recipientCentroid<-apply(scaled[!hiv[selector,'isDonor']&hiv[selector,'Fluid']=='PL'&hiv[selector,'Selection']=='UT',],2,mean)
recipientDiff<-t(t(scaled)-recipientCentroid)
recipientDist<-apply(recipientDiff^2,1,sum)

cols<-rainbow(length(unique(hiv[,'Pair ID'])),alpha=.7)
names(cols)<-sort(unique(hiv[,'Pair ID']))
pdf('out/centroidDist.pdf')
  par(mar=c(5.2,4,.1,.1))
  vpPlot(factor(hiv[selector,'fluidSelectDonor'],levels=unique(hiv$fluidSelectDonor[order(hiv$isDonor,hiv$Selection=='UT',hiv$Fluid=='PL',hiv$Selection=='A2',decreasing=TRUE)])),recipientDist,las=3,col=NA,bg=cols[as.character(hiv[selector,'Pair ID'])],pch=21,ylab='Distance to recipient samples centroid',las=2)
  legend('topright',names(cols),pch=21,pt.bg=cols,col=NA,inset=.01,title='Pair',ncol=2)
  vpPlot(factor(hiv[selector,'fluidSelectDonor'],levels=unique(hiv$fluidSelectDonor[order(hiv$isDonor,hiv$Selection=='UT',hiv$Fluid=='PL',hiv$Selection=='A2',decreasing=TRUE)])),pairRecDist,las=3,col=NA,bg=cols[as.character(hiv[selector,'Pair ID'])],pch=21,ylab='Distance to within-pair recipient centroid',las=2,cex=2)
  legend('topright',names(cols),pch=21,pt.bg=cols,col=NA,inset=.01,title='Pair',ncol=2,cex=2)
  pos<-vpPlot(factor(hiv[selector,'fluidSelectDonor'],levels=unique(hiv$fluidSelectDonor[order(hiv$isDonor,hiv$Selection=='UT',hiv$Fluid=='PL',hiv$Selection=='A2',decreasing=TRUE)])),pairRecDist,las=3,col=NA,bg=cols[as.character(hiv[selector,'Pair ID'])],pch=21,ylab='Distance to within-pair recipient centroid',las=2)
  legend('topright',names(cols),pch=21,pt.bg=cols,col=NA,inset=.01,title='Pair',ncol=2)
  text(pos,pairRecDist,hiv[selector,'Name'],cex=.2) 
dev.off()


