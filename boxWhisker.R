
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
    isLog<-targetColLog[targetCol]
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



