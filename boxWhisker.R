
if(!exists('hiv'))source('readData.R')

hiv$nameFluid<-paste(hiv$baseName,hiv$fluid)
desiredOrder<-c(
  "CH742 PL",
  "CH742 SE",
  "CH378 PL",
  "CH831 PL",
  "CH148 PL",
  "CH40 PL",
  "CH728 PL",
  "CH302 PL",
  "CH492 PL",
  "CH492 CV",
  "CH427 PL",
  "CH596 PL",
  "CH596 CV",
  "CH455 PL",
  "CH212 PL",
  "CH162 PL",
  "CH1064 PL",
  "CH848 PL"
)

for(targetCol in names(targetCols)){
  plotInfo<-do.call(rbind,lapply(desiredOrder,function(xx){
    thisDat<-hiv[hiv$nameFluid==xx&hiv$select=='UT'&!is.na(hiv[,targetCol]),targetCol]
    boxStat<-boxplot(thisDat,plot=FALSE)$stat
    if(length(thisDat)==0)thisDat<-NA
    out<-data.frame(
      'mean'=mean(thisDat),
      'geoMean'=exp(mean(log(thisDat))),
      'max'=max(thisDat),
      'min'=min(thisDat),
      'upperQuart'=boxStat[2,],
      'lowerQuart'=boxStat[4,]
    )
    return(out)
  }))
  rownames(plotInfo)<-desiredOrder
  write.csv(plotInfo,sprintf('out/boxWhisker/%s.csv',targetCol))
  pdf(sprintf('out/boxWhisker/%s.pdf',targetCol),height=5,width=5)
    par(mar=c(5.4,3,.3,.1))
    isLog<-targetColTransform[targetCol]=='log'
    ylim<-range(plotInfo,na.rm=TRUE)
    if(diff(log10(ylim))<1){
      if(log10(ylim[1])%%1<1-log10(ylim[2])%%1)ylim[1]<-10^floor(log10(ylim[1]))
      else ylim[2]<-10^ceiling(log10(ylim[2]))
    }
    plot(1,1,type='n',xlab='',ylab=targetCols[targetCol],xlim=c(1,nrow(plotInfo)),ylim=ylim,xaxt='n',mgp=c(2,1,0),las=1,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'))
    if(isLog)logAxis(2,las=1,mgp=c(3,.5,0))
    segments(1:nrow(plotInfo),plotInfo$max,1:nrow(plotInfo),plotInfo$min)
    rect(1:nrow(plotInfo)-.2,plotInfo$upperQuart,1:nrow(plotInfo)+.2,plotInfo$lowerQuart)
    means<-plotInfo[,ifelse(isLog,'geoMean','mean')]
    segments(1:nrow(plotInfo)-.2,means,1:nrow(plotInfo)+.2,means)
    axis(1,1:nrow(plotInfo),rownames(plotInfo),las=3)
  dev.off()
}


cols<-rainbow.lab(length(unique(hiv$Pair.ID)),alpha=.2)
cols2<-rainbow.lab(length(unique(hiv$Pair.ID)),alpha=.7)
names(cols)<-names(cols2)<-unique(hiv$Pair.ID)
for(targetCol in names(targetCols)){
  selector<-!is.na(hiv[,targetCol])
  pos<-seq(-.3,.3,length.out=max(hiv$Pair.ID))
  catPos<-structure(1:length(unique(hiv$fluidSelectDonor)),names=unique(hiv$fluidSelectDonor[order(hiv$donor,hiv$fluid=='PL',hiv$select=='UT',hiv$select=='A2',decreasing=TRUE)]))
  pdf(sprintf('out/boxWhisker/7line_%s.pdf',targetCol),width=8)
    par(mar=c(5,5,.1,.1))
    plot(1,1,type='n',xlim=c(.5,length(unique(hiv$fluidSelectDonor))),ylim=range(hiv[selector,targetCol]),ylab=targetCols[targetCol],log='y',las=1,xaxt='n',yaxt='n',xlab='',mgp=c(4,1,0))
    logAxis(2,las=1,addExtra=TRUE,mgp=c(3,.7,0))
    axis(1,catPos,names(catPos),las=2)
    rect(seq(1.5,max(catPos),2),10^par('usr')[3],seq(2.5,max(catPos)+.5,2),10^par('usr')[4],col='#00000011',border=NA)
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],hiv[selector,'fluidSelectDonor'],range))
    xPos<-catPos[hiv[selector,'fluidSelectDonor']]+pos[hiv[selector,'Pair.ID']]
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],xPos,range))
    segments(as.numeric(rownames(ranges)),ranges[,1],as.numeric(rownames(ranges)),ranges[,2],col='#00000033')
    points(xPos+offsetX(log10(hiv[selector,targetCol]),xPos,width=.04),hiv[selector,targetCol],cex=.5,pch=21,col=NA,bg=cols2[hiv[selector,'Pair.ID']])
  dev.off()
}

exp(mean(log(hiv[hiv$baseName=='CH492'&hiv$select=='A2'&hiv$fluid=='PL','IFNa2.Pooled.Donor.cells.IC50..pg..ml'])))/
exp(mean(log(hiv[hiv$baseName=='CH492'&hiv$select=='UT'&hiv$fluid=='PL','IFNa2.Pooled.Donor.cells.IC50..pg..ml'])))

exp(mean(log(hiv[hiv$baseName=='CH492'&hiv$select=='BE'&hiv$fluid=='PL','IFNbeta.Pooled.Donor.cells.IC50..pg.ml'])))/
exp(mean(log(hiv[hiv$baseName=='CH492'&hiv$select=='UT'&hiv$fluid=='PL','IFNbeta.Pooled.Donor.cells.IC50..pg.ml'])))
