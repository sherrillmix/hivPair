library(levenR)
library(vipor)
library(dnar)
library(dnaplotr)
library(glmnet)
library(glmnetPlotR)
library(parallel)
#library(beeswarm)

hiv<-read.csv('data/data.csv',stringsAsFactors=FALSE)
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

align<-read.fa('data/AlignSeq.nt.fasta.gz')
seqMat<-seqSplit(align$seq)[-1,]
refPos<-cumsum(seqMat[1,]!='-')
colnames(seqMat)<-sprintf('pos%s.%s',refPos,ave(refPos,refPos,FUN=function(x)1:length(x)))
nBase<-apply(seqMat,2,function(x)sum(x!='-'))
startEnd<-range(which(nBase/nrow(seqMat)>.9))
seqMat<-seqMat[,nBase>10&1:ncol(seqMat) %in% startEnd[1]:startEnd[2]]
onlyDiff<-seqMat[,apply(seqMat,2,function(x)length(unique(x[x!='-']))>1)]
png('out/diff.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotDNA(apply(onlyDiff,1,paste,collapse='')[hiv$donor],groups=paste(hiv$Pair.ID..,hiv$select)[hiv$donor])
dev.off()
png('out/align.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotDNA(apply(seqMat,1,paste,collapse=''),groups=paste(hiv$Pair.ID..,hiv$donor))
dev.off()



