library(vipor)
if(!exists('hiv'))source('readData.R')

selector<-!apply(is.na(hiv[,names(goodTargetCols)]),1,any)
tmp<-hiv[selector,names(goodTargetCols)]
rownames(tmp)<-hiv[selector,'Renamed']
pca<-prcomp(tmp,scale.=TRUE)
pcaPoints<-t(t(pca$x)/pca$sdev/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
importance<-summary(pca)$importance[2,]
cols<-rainbow.lab(length(unique(hiv$fluidSelectDonor)),alpha=.6)
#pch<-structure(c(21,22,22),names=c('PL','SE','CV'))
pch<-structure(c(21,21,21),names=c('PL','SE','CV'))
cols2<-rainbow.lab(length(unique(hiv$fluidSelectDonor)),alpha=.8)
names(cols)<-names(cols2)<-sort(unique(hiv$fluidSelectDonor))
pdf('out/pca.pdf')
  for(select in list(1:2,2:3,3:4,4:5)){
    plot(
      1,1,type='n',las=1,
      xlim=range(-pcaPoints[,select[1]]),ylim=range(-pcaPoints[,select[2]]),
      xlab=sprintf('Principal component %d (%d%% of variance)',select[1],round(importance[select[1]]*100)),
      ylab=sprintf('Principal component %d (%d%% of variance)',select[2],round(importance[select[2]]*100))
    )
    points(-pcaPoints[,select[1]],-pcaPoints[,select[2]],bg=cols[as.character(hiv[selector,'fluidSelectDonor'])],col=cols2[as.character(hiv[selector,'fluidSelectDonor'])],pch=21,cex=1.5)
    legend('topright',names(cols),pch=21,col=cols2,pt.bg=cols,inset=.01,pt.cex=1.5)
    biplot(pca,choices=select,cex=.25)
  }
dev.off()

scaled<-apply(hiv[selector,names(goodTargetCols)],2,function(x)(x-mean(x))/sd(x))
sampleNames<-withAs(xx=hiv[selector,],paste(xx$sampleFluidSelect,ave(xx$sampleFluidSelect,xx$sampleFluidSelect,FUN=function(x)1:length(x))))
rownames(scaled)<-sampleNames
dists<-dist(scaled)
distMat<-as.matrix(dists)


ids<-withAs(xx=hiv[selector,],list(
  donor=sampleNames[xx$donor&xx$select=='UT'&xx$fluid=='PL']
  ,donorGenital=sampleNames[xx$donor&xx$select=='UT'&xx$fluid!='PL']
  ,donorAlpha=sampleNames[xx$donor&xx$select=='A2'&xx$fluid=='PL']
  ,donorBeta=sampleNames[xx$donor&xx$select=='BE'&xx$fluid=='PL']
  ,recAlpha=sampleNames[!xx$donor&xx$select=='A2'&xx$fluid=='PL']
  ,recBeta=sampleNames[!xx$donor&xx$select=='BE'&xx$fluid=='PL']
  ,recipient=sampleNames[!xx$donor&xx$select=='UT']
))

distList<-lapply(ids,function(xx,distMat,recipientIds){
  select<-distMat[xx,recipientIds]
  dists<-select[upper.tri(select)]
  return(dists)
},distMat,ids[['recipient']])
distList<-distList[order(sapply(distList,mean),decreasing=TRUE)]

pdf('out/dists.pdf')
par(mar=c(7.2,4,.1,.1))
  vpPlot(factor(rep(names(distList),sapply(distList,length)),levels=names(distList)),unlist(distList),las=3,cex=.5,col=NA,bg='#00000066',pch=21,ylab='Distance to recipient samples')
dev.off()

pairRecDist<-unlist(ave(split(cbind(scaled,withAs(xx=hiv[selector,],!xx$donor&xx$fluid=='PL'&xx$select=='UT')),1:nrow(scaled)),hiv[selector,'Pair.ID'],FUN=function(df){
  df<-do.call(rbind,df)
  select<-df[,ncol(df)]==1
  df<-df[,-ncol(df)]
  recipientCentroid<-apply(df[select,],2,mean)
  recipientDiff<-t(t(df)-recipientCentroid)
  return(apply(recipientDiff^2,1,sum))
}))

recipientCentroid<-apply(scaled[withAs(xx=hiv[selector,],!xx$donor&xx$fluid=='PL'&xx$select=='UT'),],2,mean)
recipientDiff<-t(t(scaled)-recipientCentroid)
recipientDist<-apply(recipientDiff^2,1,sum)

cols<-rainbow.lab(length(unique(hiv$Pair.ID)),alpha=.7)
names(cols)<-sort(unique(hiv$Pair.ID))
pdf('out/centroidDist.pdf')
  par(mar=c(5.2,4,.1,.1))
  vpPlot(factor(hiv[selector,'fluidSelectDonor'],levels=unique(hiv$fluidSelectDonor[order(hiv$donor,hiv$select=='UT',hiv$fluid=='PL',hiv$select=='A2',decreasing=TRUE)])),recipientDist,las=3,col=NA,bg=cols[as.character(hiv[selector,'Pair.ID'])],pch=21,ylab='Distance to recipient samples centroid',las=2)
  legend('topright',names(cols),pch=21,pt.bg=cols,col=NA,inset=.01,title='Pair',ncol=2)
  vpPlot(factor(hiv[selector,'fluidSelectDonor'],levels=unique(hiv$fluidSelectDonor[order(hiv$donor,hiv$select=='UT',hiv$fluid=='PL',hiv$select=='A2',decreasing=TRUE)])),pairRecDist,las=3,col=NA,bg=cols[as.character(hiv[selector,'Pair.ID'])],pch=21,ylab='Distance to within-pair recipient centroid',las=2)
  legend('topright',names(cols),pch=21,pt.bg=cols,col=NA,inset=.01,title='Pair',ncol=2)
dev.off()

pdf('out/tree.pdf',height=50,width=20);par(mar=c(4,1,1,10));plot(as.dendrogram(hclust(dists)),horiz=TRUE);dev.off()

