
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

targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Pooled donor\nReplicative capacity (day 7 p24)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml.'='IFNbeta IC50 (pg/ml)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'='IFNa2 IC50 (pg/ml)'
)

for(targetCol in names(targetCols)){
  plotInfo<-do.call(rbind,lapply(desiredOrder,function(xx){
    thisDat<-hiv[hiv$nameFluid==xx&hiv$select=='UT'&!is.na(hiv[,targetCol]),targetCol]
    boxStat<-boxplot(thisDat,plot=FALSE)$stat
    out<-data.frame(
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
    plot(1,1,type='n',xlab='',ylab=targetCols[targetCol],xlim=c(1,nrow(plotInfo)),ylim=range(plotInfo),xaxt='n',mgp=c(2,1,0),las=1,log='y',yaxt='n')
    logAxis(2,las=1)
    segments(1:nrow(plotInfo),plotInfo$max,1:nrow(plotInfo),plotInfo$min)
    rect(1:nrow(plotInfo)-.2,plotInfo$upperQuart,1:nrow(plotInfo)+.2,plotInfo$lowerQuart)
    segments(1:nrow(plotInfo)-.2,plotInfo$geoMean,1:nrow(plotInfo)+.2,plotInfo$geoMean)
    axis(1,1:nrow(plotInfo),rownames(plotInfo),las=3)
  dev.off()
}



