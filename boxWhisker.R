library(vipor)

if(!exists('hiv'))source('readData.R')
if(!dir.exists(file.path('out','boxWhisker')))dir.create(file.path('out','boxWhisker'))

hiv$nameFluid<-paste(hiv$baseName,hiv$Fluid)
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

#generate box whisker plots
for(targetCol in targetCols){
  #calculate stats
  plotInfo<-do.call(rbind,lapply(desiredOrder,function(xx){
    thisDat<-hiv[hiv$nameFluid==xx&hiv$Selection=='UT'&!is.na(hiv[,targetCol]),targetCol]
    #making sure quantile agrees with boxplot quantile calculation
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
  pdf(file.path('out','boxWhisker',sprintf('%s.pdf',gsub('/','_',targetCol))),height=5,width=5)
    par(mar=c(5.4,4,.3,.1))
    isLog<-targetColTransform[targetCol]=='log'
    ylim<-range(plotInfo,na.rm=TRUE)
    if(diff(log10(ylim))<1){
      if(log10(ylim[1])%%1<1-log10(ylim[2])%%1)ylim[1]<-10^floor(log10(ylim[1]))
      else ylim[2]<-10^ceiling(log10(ylim[2]))
    }
    plot(1,1,type='n',xlab='',ylab=targetCols[targetCol],xlim=c(1,nrow(plotInfo)),ylim=ylim,xaxt='n',mgp=c(2,1,0),las=1,log=ifelse(isLog,'y',''))
    #add whiskers
    segments(1:nrow(plotInfo),plotInfo$max,1:nrow(plotInfo),plotInfo$min)
    #add box
    rect(1:nrow(plotInfo)-.2,plotInfo$upperQuart,1:nrow(plotInfo)+.2,plotInfo$lowerQuart)
    means<-plotInfo[,ifelse(isLog,'geoMean','mean')]
    #mean segment
    segments(1:nrow(plotInfo)-.2,means,1:nrow(plotInfo)+.2,means)
    axis(1,1:nrow(plotInfo),rownames(plotInfo),las=3)
  dev.off()
}


#add transparency to colors
cols<-sprintf('%s33',pairColors)
cols2<-sprintf('%sB3',pairColors)

#generate seven line plots
for(targetCol in targetCols){
  selector<-!is.na(hiv[,targetCol])
  #quick function to calculate out position for each line of points
  pos<-as.numeric(ave(apply(hiv[selector,c('Donor/Recipient','Pair ID')],1,paste,collapse=''),hiv$fluidSelectDonor[selector],FUN=function(pair){
    nPos<-length(unique(pair))
    step<-.13
    if(nPos<4)step<-.2
    if(nPos>7)step<-.11
    xStart<--(nPos-1)*step/2
    pairPos<-structure(xStart+seq(0,nPos*step,step),names=sort(unique(pair)))
    pairPos[as.character(pair)]
  }))
  catPos<-structure(1:length(unique(hiv$fluidSelectDonor)),names=unique(hiv$fluidSelectDonor[order(hiv$isDonor,hiv$Selection=='UT',!hiv$isGenital,hiv$Selection=='A2',decreasing=TRUE)]))
  pdf(file.path('out','boxWhisker',sprintf('7line_%s.pdf',sub('/','_',targetCol))),width=8,height=4)
    par(mar=c(3.2,5,.1,.1))
    logY<-ifelse(targetColTransform[targetCol]=='log','y','')
    plot(1,1,type='n',xlim=c(.4,length(unique(hiv$fluidSelectDonor))+.6),ylim=range(hiv[selector,targetCol]),ylab=targetCol,log=logY,las=1,xaxt='n',xlab='',mgp=c(4,1,0),xaxs='i')
    par(lheight=.8)
    axis(1,catPos,gsub(' ','\n',names(catPos)),las=1,padj=1,mgp=c(3,0,0))
    if(logY=='y')rect(seq(1.5,max(catPos),2),10^par('usr')[3],seq(2.5,max(catPos)+.5,2),10^par('usr')[4],col='#00000011',border=NA)
    else rect(seq(1.5,max(catPos),2),par('usr')[3],seq(2.5,max(catPos)+.5,2),par('usr')[4],col='#00000011',border=NA)
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],hiv[selector,'fluidSelectDonor'],range))
    xPos<-catPos[hiv[selector,'fluidSelectDonor']]+pos
    ranges<-do.call(rbind,tapply(hiv[selector,targetCol],xPos,range))
    #add segment for each patient
    segments(as.numeric(rownames(ranges)),ranges[,1],as.numeric(rownames(ranges)),ranges[,2],col='#00000011')
    #plot point for each virus offset slightly using vipor library
    points(xPos+offsetX(log10(hiv[selector,targetCol]),xPos,width=.05),hiv[selector,targetCol],cex=1.3,pch=21,col=NA,bg=cols2[hiv[selector,'Pair ID']])
  dev.off()
}

