library(levenR)
library(vipor)
library(dnar)
#library(beeswarm)

hiv<-read.csv('data.csv',stringsAsFactors=FALSE)
hiv<-hiv[,apply(hiv,2,function(x)!all(is.na(x)))]
hiv$seq<-toupper(gsub('[ \n]','',hiv$Sequence))
hiv$select<-sapply(strsplit(hiv$Renamed,'\\.'),'[',4)
hiv$donorRec<-sub(' ','-',hiv$Donor.or.Recipient)
hiv$sample<-paste(hiv$donorRec,hiv$Pair.ID..)
hiv$donor<-grepl('Donor',hiv$Donor.or.Recipient)
hiv$sampleSelect<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$select)

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
  selector<-!is.na(thisData)&hiv$select=='UT'
  par(mar=c(6.5,4,.1,.1),mgp=c(3,.7,0))
  plot(1,1,las=2,ylab=colnames(hiv)[ii],xlim=range(labPos),ylim=range(thisData,na.rm=TRUE),xaxt='n',xlab='',type='n')
  axis(1,labPos,uniqSample,las=2)
  xOffset<-offsetX(thisData[selector],labs[selector])
  rect(c(1,labPos[changes])-.75,par('usr')[3],c(labPos[changes],length(pairId)+1)-.75,par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=.9,col=NA,bg=cols[hiv$donor[selector]+1],pch=21)
}
dev.off()

uniqSample<-unique(hiv$sampleSelect[order(hiv$Pair.ID..,hiv$donorRec,c('UT'=1,'A2'=2,'BE'=3)[hiv$select])])
pairId<-sub('^[A-Za-z0-9-]+ ([0-9]) [A-Z0-9]+','\\1',uniqSample)
labs<-factor(hiv$sampleSelect,levels=uniqSample)
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

