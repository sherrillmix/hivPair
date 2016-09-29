library(vipor)
if(!exists('hiv'))source('readData.R')

selector<-!apply(is.na(hiv[,names(targetCols)]),1,any)
pca<-prcomp(hiv[selector,names(targetCols)],scale.=TRUE)
pcaPoints<-t(t(pca$x)/pca$sdev/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
importance<-summary(pca)$importance[2,]

cols<-rainbow.lab(length(unique(hiv$donor)),alpha=.8)
names(cols)<-unique(hiv$donor)
pch<-structure(c(21,22,22),names=c('PL','SE','CV'))
col=ifelse(hiv[selector,'select']=='BE','black',ifelse(hiv[selector,'select']=='A2','red','#00000055'))

pdf('out/pca.pdf')
  for(select in list(1:2,2:3,3:4,4:5)){
    plot(
      1,1,type='n',las=1,
      xlim=range(pcaPoints[,select[1]]),ylim=range(pcaPoints[,select[2]]),
      xlab=sprintf('Principal component %d (%d%% of variance)',select[1],round(importance[select[1]]*100)),
      ylab=sprintf('Principal component %d (%d%% of variance)',select[2],round(importance[select[2]]*100))
    )
    for(ii in c('UT','A2','BE')){
      selectSelector<-hiv[selector,'select']==ii
      points(pcaPoints[selectSelector,],bg=cols[as.character(hiv[selector,'donor'])][selectSelector],pch=pch[hiv[selector,'fluid']][selectSelector],col=col[selectSelector])
    }
    legend('bottomleft',c('Genital','Plasma','IFNa2 select','IFNb select','Donor','Recipient'),pch=c(22,21,rep(21,4)),col=c('#00000055','#00000055','red','black','#00000055','#00000055'),pt.bg=rep(cols[c('TRUE','FALSE')],c(5,1)),inset=.01)
  }
  #plot(pca)
dev.off()

scaled<-apply(hiv[selector,names(targetCols)],2,function(x)(x-mean(x))/sd(x))
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

pdf('out/tree.pdf',height=50,width=20);par(mar=c(4,1,1,10));plot(as.dendrogram(hclust(dists)),horiz=TRUE);dev.off()

