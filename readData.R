library(levenR)
library(vipor)
library(dnar)
#library(beeswarm)

hiv<-read.csv('data.csv',stringsAsFactors=FALSE)
hiv<-hiv[,apply(hiv,2,function(x)!all(is.na(x)))]
hiv$seq<-toupper(gsub('[ \n]','',hiv$Sequence))
hiv$select<-sapply(strsplit(hiv$Renamed,'\\.'),'[',4)
hiv$fluid<-sapply(strsplit(hiv$Renamed,'\\.'),'[',2)
hiv$donorRec<-sub(' ','-',hiv$Donor.or.Recipient)
hiv$sample<-paste(hiv$donorRec,hiv$Pair.ID..)
hiv$donor<-grepl('Donor',hiv$Donor.or.Recipient)
hiv$sampleSelect<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$select)
hiv$sampleFluid<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$fluid)
write.fa(hiv$name,hiv$seq,'hiv.fa')

varCols<-which(colnames(hiv)=='Replicative.capacity.Single.Donor.p24.d7'):which(colnames(hiv)=='Bnaber.IC50')
pdf('out/pairs.pdf',width=12,height=12)
thisData<-hiv[,varCols]
colnames(thisData)<-gsub('\\.+','\n',colnames(thisData))
plot(thisData,col='#00000077',pch='.',cex=2.5)
dev.off()

#align<-levenAlign(hiv$seq,hiv$seq[1])



uniqSample<-unique(hiv$sample[order(hiv$Pair.ID..,hiv$Donor.or.Recipient)])
pairId<-sub('.* ','',uniqSample)
labs<-factor(hiv$sample,levels=uniqSample)
cols<-c('#FF000066','#0000FF66')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
changes<-which(isChange)
labPos<-1:length(uniqSample)+cumsum(isChange)*.5
xPos<-labPos[as.numeric(labs)]
pdf('out/untreated_vpPlot.pdf',height=4,width=8)
for(ii in varCols){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)&hiv$select=='UT'&hiv$fluid=='PL'
  par(mar=c(6.5,4,.1,.1),mgp=c(3,.7,0))
  plot(1,1,las=2,ylab=colnames(hiv)[ii],xlim=range(labPos),ylim=range(thisData,na.rm=TRUE),xaxt='n',xlab='',type='n')
  axis(1,labPos,uniqSample,las=2)
  xOffset<-offsetX(thisData[selector],labs[selector])
  rect(c(1,labPos[changes])-.75,par('usr')[3],c(labPos[changes],length(pairId)+1)-.75,par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=.9,col=NA,bg=cols[hiv$donor[selector]+1],pch=21)
}
dev.off()






uniqSample<-unique(hiv$sampleSelect[order(hiv$Pair.ID..,hiv$donorRec,c('UT'=1,'A2'=2,'BE'=3)[hiv$select])])
pairId<-sub('^[A-Za-z0-9-]+ ([0-9]) [A-Z0-9]+','\\1',uniqSample) labs<-factor(hiv$sampleSelect,levels=uniqSample)
cols<-c('#FF000066','#0000FF66')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
donorRec<-sub(' .+','',uniqSample)
treat<-sub('.* ','',uniqSample)
isChange2<-c(FALSE,donorRec[-1]!=donorRec[-length(pairId)])
changes<-which(isChange)
changes2<-which(c(FALSE,treat[-1]!=treat[-length(pairId)]))
labPos<-1:length(uniqSample)+cumsum(isChange)*.5+cumsum(isChange2*.5)
xPos<-labPos[as.numeric(labs)]
cols<-rainbow.lab(10,alpha=.7)[c(1:3,8:10)]
names(cols)<-unique(paste(hiv$donor,hiv$select))

pdf('out/selects_vpPlot.pdf',height=4,width=8)
for(ii in varCols){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)
  par(mar=c(7.3,4,.3,.1),mgp=c(3,.7,0))
  plot(1,1,las=2,ylab=colnames(hiv)[ii],xlim=range(labPos)+c(-.7,.7),ylim=range(thisData,na.rm=TRUE),xaxt='n',xlab='',type='n',log='y',xaxs='i')
  axis(1,labPos,uniqSample,las=2)
  xOffset<-offsetX(log10(thisData[selector]),labs[selector],width=.4)
  #xOffset<-ave(thisData[selector],labs[selector],FUN=function(x)swarmx(rep(0,length(x)),x,log='y',cex=.9)$x)
  rect(c(1,labPos[changes])-1,10^par('usr')[3],c(labPos[changes],length(pairId)+1)-1,10^par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  #rect(c(1,labPos[changes2])-.5,10^par('usr')[3],c(labPos[changes2],length(pairId)+1)-.5,10^par('usr')[4],col=rep(c(NA,'#00000006'),length.out=length(changes2)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=.9,col=NA,bg=cols[paste(hiv$donor,hiv$select)[selector]],pch=21)
}
dev.off()


#Env  content ( Env :RT ratios)
#Infectivity
#Replicative capacity
#IFN alpha IC50
selectVars<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pool.Donor.p24.d7'='Pooled donor\nReplicative capactity (day 7 p24)',
  'IFNa2.PD.IC50..U.ml.'='Pooled donor\nIFNa2 IC50 (U/ml)'
)
  #'Replicative.capacity.Single.Donor.p24.d7'='Single donor\nReplicative capactity\n(day 7 p24)',
  #'IFNa2.SD.IC50..U.ml.'='Single donor\nIFNa2 IC50 (U/ml)',

logAxis<-function(x,axisNum=2,spreadRange=1.3){
  minX<-min(log10(x),na.rm=TRUE) 
  maxX<-max(log10(x),na.rm=TRUE) 
  allTicks<-unlist(lapply(floor(minX):ceiling(maxX),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<10^maxX*spreadRange & allTicks>10^minX/spreadRange]
  axis(axisNum,allTicks,rep('',length(allTicks)),tcl=-.2)
  prettyY<-seq(ceiling(log10(min(x,na.rm=TRUE))),floor(log10(max(x,na.rm=TRUE))),1)
  if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(5),prettyY-log10(2)))
  if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(2),prettyY-log10(5)))
  axis(axisNum,10^prettyY,10^prettyY,las=2)
}
uniqSample<-unique(hiv$sampleFluid[order(hiv$Pair.ID..,hiv$Donor.or.Recipient,c('PL'=1,'CV'=2,'SE'=3)[hiv$fluid])])
pairId<-sub('.* ([0-9]+) [A-Z]+','\\1',uniqSample)
labs<-factor(hiv$sampleFluid,levels=uniqSample)
cols<-c('TRUE PL'='#FF000055','TRUE CV'='#0000FF55','TRUE SE'='#0000FF55','FALSE PL'='#FF770055')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
changes<-which(isChange)
labPos<-1:length(uniqSample)+cumsum(isChange)*.5
xPos<-labPos[as.numeric(labs)]
pdf('out/selectVar.pdf',height=8,width=7)
layout(matrix(c(length(selectVars)+2,1:(length(selectVars)+1)),ncol=1),height=c(.05,rep(1,length(selectVars)),.6))
for(ii in names(selectVars)){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)&hiv$select=='UT'
  par(mar=c(0,4.8,0,.2),mgp=c(3,.7,0),lheight=.8)
  plot(1,1,las=2,ylab=selectVars[ii],xlim=range(labPos),ylim=range(thisData,na.rm=TRUE)*c(.8,1.2),xaxt='n',xlab='',type='n',log='y',yaxt='n')
  logAxis(thisData)
  xOffset<-offsetX(thisData[selector],labs[selector],width=.45,varwidth=TRUE)
  rect(c(1,labPos[changes])-.75,10^par('usr')[3],c(labPos[changes],length(pairId)+1)-.75,10^par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=1,col=NA,bg=cols[paste(hiv$donor,hiv$fluid)[selector]],pch=21)
  #axis(1,labPos,rep('',length(labPos)))
}
axis(1,labPos,uniqSample,las=2)
dev.off()



thisData<-hiv[hiv$select=='UT'&hiv$fluid=='PL',]
infectData<-tapply(thisData$Infectivity.RLU.pg.RT...T1249.,list(thisData$Donor.or.Recipient,thisData$Pair.ID..),median,na.rm=TRUE)
infectData25<-tapply(thisData$Infectivity.RLU.pg.RT...T1249.,list(thisData$Donor.or.Recipient,thisData$Pair.ID..),quantile,probs=c(.25),na.rm=TRUE)
infectData75<-tapply(thisData$Infectivity.RLU.pg.RT...T1249.,list(thisData$Donor.or.Recipient,thisData$Pair.ID..),quantile,probs=c(.75),na.rm=TRUE)
ylim<-range(infectData,na.rm=TRUE)
xpos<-((1:ncol(infectData))/ncol(infectData)-.5)*.2
#xpos<-order(ifelse(is.na(infectData[2,]),infectData[3,],infectData[2,]))
cols<-rainbow.lab(ncol(infectData),alpha=.7)
pdf('out/pairInfect.pdf')
  plot(1,1,type='n',xlim=c(.75,2.25),ylim=ylim,ylab='Infectivity (RLU/pg RT)',xlab='',xaxt='n',log='y',las=1)
  axis(1,1:2,c('Donor','Recipient'))
  for(ii in 1:ncol(infectData)){
    segments(1+xpos[ii],infectData['Donor',ii],2+xpos[ii],infectData[rownames(infectData)!='Donor',ii],col=cols[ii],lty=2)
    segments(1+xpos[ii],infectData25['Donor',ii],1+xpos[ii],infectData75['Donor',ii],col=cols[ii])
    segments(2+xpos[ii],infectData25[rownames(infectData)!='Donor',ii],2+xpos[ii],infectData75[rownames(infectData)!='Donor',ii],col=cols[ii])
  }
dev.off()

