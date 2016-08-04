library(levenR)
library(vipor)

hiv<-read.csv('data.csv',stringsAsFactors=FALSE)
hiv<-hiv[,apply(hiv,2,function(x)!all(is.na(x)))]
hiv$seq<-toupper(gsub('[ \n]','',hiv$Sequence))
hiv$select<-sapply(strsplit(hiv$Renamed,'\\.'),'[',4)
hiv$sample<-paste(hiv$Donor.or.Recipient,hiv$Pair.ID..)
hiv$donor<-grepl('Donor',hiv$Donor.or.Recipient)

varCols<-which(colnames(hiv)=='Replicative.capacity.Single.Donor.p24.d7'):which(colnames(hiv)=='Bnaber.IC50')
pdf('out/pairs.pdf',width=12,height=12)
thisData<-hiv[,varCols]
colnames(thisData)<-gsub('\\.+','\n',colnames(thisData))
plot(thisData,col='#00000077',pch='.',cex=2.5)
dev.off()

align<-levenAlign(hiv$seq,hiv$seq[1])



uniqSample<-unique(hiv$sample[order(hiv$Pair.ID..,hiv$Donor.or.Recipient)])
pairId<-sub('.* ','',uniqSample)
labs<-factor(hiv$sample,levels=uniqSample)
cols<-c('#FF000099','#0000FF99')
changes<-which(c(pairId[-1]!=pairId[-length(pairId)]))
pdf('out/vpPlot.pdf',height=6,width=8)
for(ii in varCols){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)
  par(mar=c(6.5,4,.1,.1),mgp=c(3,.7,0))
  vpPlot(labs[selector],thisData[selector],pch='.',cex=3,las=2,ylab=colnames(hiv)[ii],col=cols[hiv$donor[selector]+1])
  rect(c(0,changes)+.5,par('usr')[3],c(changes,length(pairId))+.5,par('usr')[4],col=rep(c(NA,'#00000033'),length.out=length(changes)),border=NA)
}
dev.off()

