
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
    print(isLog)
    plot(1,1,type='n',xlab='',ylab=targetCols[targetCol],xlim=c(1,nrow(plotInfo)),ylim=range(plotInfo),xaxt='n',mgp=c(2,1,0),las=1,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'))
    if(isLog)logAxis(2,las=1)
    segments(1:nrow(plotInfo),plotInfo$max,1:nrow(plotInfo),plotInfo$min)
    rect(1:nrow(plotInfo)-.2,plotInfo$upperQuart,1:nrow(plotInfo)+.2,plotInfo$lowerQuart)
    means<-plotInfo[,ifelse(isLog,'geoMean','mean')]
    segments(1:nrow(plotInfo)-.2,means,1:nrow(plotInfo)+.2,means)
    axis(1,1:nrow(plotInfo),rownames(plotInfo),las=3)
  dev.off()
}


for(targetCol in names(targetCols)){
  selector<-!is.na(hiv[,targetCol])
  pos<-seq(-.2,.2,length.out=max(hiv$Pair.ID..))
  catPos<-structure(1:length(unique(hiv$fluidSelectDonor)),names=unique(hiv$fluidSelectDonor[order(hiv$donor,hiv$fluid,hiv$select=='PL',decreasing=TRUE)]))
  pdf(sprintf('out/boxWhisker/7line_%s.pdf',targetCol))
    plot(1,1,type='n',xlim=c(.5,length(unique(hiv$fluidSelectDonor))),ylim=range(hiv[selector,targetCol]),ylab=targetCols[targetCol],log='y',las=1,xaxt='n',yaxt='n',xlab='')
    logAxis(2)
    axis(1,catPos,names(catPos),las=2)
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],hiv[selector,'fluidSelectDonor'],range))
    xPos<-catPos[hiv[selector,'fluidSelectDonor']]+pos[hiv[selector,'Pair.ID..']]
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],xPos,range))
    segments(as.numeric(rownames(ranges)),ranges[,1],as.numeric(rownames(ranges)),ranges[,2],col='#00000044')
    points(xPos,hiv[selector,targetCol],cex=.5,pch=21,col=NA,bg='#00000088')
  dev.off()
}


