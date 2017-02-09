if(!exists('hiv'))source('readData.R')
iso<-read.csv('data/IsolationSummary.csv',stringsAsFactors=FALSE)
iso$donor<-sapply(iso$ID,function(x)hiv[hiv$baseName==x,'donor'][1])
iso$pair<-sapply(iso$ID,function(x)hiv[hiv$baseName==x,'Pair.ID'][1])
iso$infectPerVirus<-iso$X.positiove.wells..d20/(iso$vRNA.well*iso$X..24.well.plates*24)
iso$select<-sub(' +$','',iso$Selection)

iso[order(iso$pair,iso$donor),c('donor','pair','select','infectPerVirus')]

isoTable<-tapply(iso$infectPerVirus,list(paste(iso$pair,ifelse(iso$donor,'D','R'),sep=''),iso$select),mean)
isoTable<-isoTable[,order(colnames(isoTable)!='UT',colnames(isoTable)!='A2')]
cols<-rainbow.lab(3)
names(cols)<-unique(iso$select)
prettyNames<-c('UT'="Untreated",'A2'="IFNa2 treated","BE"="IFNb treated")
pdf('out/iso.pdf',width=6,height=4)
  par(mar=c(2,3,.1,.1))
  barplot(t(isoTable*100),beside=TRUE,ylab="Positive wells/viral genome (%)",col=cols[colnames(isoTable)],las=1,mgp=c(2,.8,0))
  legend('topright',prettyNames[names(cols)],fill=cols)
dev.off()
