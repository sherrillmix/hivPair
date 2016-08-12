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
hiv$seqId<-sprintf('seq%d',1:nrow(hiv))
hiv$baseName<-sub('\\..*$','',hiv$Renamed)
donors<-with(hiv[hiv$donor,],tapply(baseName,Pair.ID..,unique))
recs<-with(hiv[!hiv$donor,],tapply(baseName,Pair.ID..,unique))
pairNames<-mapply(function(x,y)paste(paste(x,collapse='/'),paste(y,collapse='/'),sep='-'),donors,recs)
write.fa(hiv$seqId,hiv$seq,'hiv.fa')

refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
align<-read.fa('data/AlignSeq.nt.fasta.gz')
rownames(align)<-align$name
seqMat<-seqSplit(align$seq)
rownames(seqMat)<-align$name
refPos<-cumsum(seqMat[refName,]!='-')
seqMat<-seqMat[hiv$seqId,]
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
plotDNA(align[hiv$seqId,'seq'],groups=paste(hiv$Pair.ID..,hiv$donor))
dev.off()


aa<-readFaDir('data','^[^LA].*aa.fasta.gz$') #don't read LTR since no AA or AlignSeq
bigAA<-tapply(aa$seq,aa$name,paste,collapse='')
aaMat<-seqSplit(bigAA)
rownames(aaMat)<-names(bigAA)
refPos<-cumsum(aaMat[refName,]!='-')
aaMat<-aaMat[hiv$seqId,]
colnames(aaMat)<-sprintf('pos%s.%s',refPos,ave(refPos,refPos,FUN=function(x)1:length(x)))
nBase<-apply(aaMat,2,function(x)sum(x!='-'))
seqMat<-seqMat[,nBase>10]
onlyDiffAA<-aaMat[,apply(aaMat,2,function(x)length(unique(x[x!='-']))>1)]


png('out/diffAA.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotAA(apply(onlyDiffAA,1,paste,collapse='')[hiv$donor],groups=paste(hiv$Pair.ID..,hiv$select)[hiv$donor])
dev.off()
png('out/alignAA.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotAA(apply(aaMat,1,paste,collapse=''),groups=paste(hiv$Pair.ID..,hiv$donor))
dev.off()


