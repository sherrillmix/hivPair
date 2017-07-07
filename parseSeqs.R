if(!exists('hiv'))source("readData.R")
refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
align<-read.fa('data/AlignSeq.nt.fasta.gz')
rownames(align)<-align$name
rawSeqMat<-seqSplit(align$seq)
rownames(rawSeqMat)<-align$name
refPos<-cumsum(rawSeqMat[refName,]!='-')
rawSeqMat<-rawSeqMat[hiv$seqId,]
if(any(degap(apply(rawSeqMat,1,paste,collapse=''))!=hiv$seq))stop(simpleError('hiv metadata and alignment mismatched'))
colnames(rawSeqMat)<-sprintf('nt%s.%s',refPos,ave(refPos,refPos,FUN=function(x)1:length(x)))
nBase<-apply(rawSeqMat,2,function(x)sum(x!='-'))
startEnd<-range(which(nBase/nrow(rawSeqMat)>.9))
seqMat<-rawSeqMat[,nBase>10&1:ncol(rawSeqMat) %in% startEnd[1]:startEnd[2]]
onlyDiff<-seqMat[,apply(seqMat,2,function(x)length(unique(x[x!='-']))>1)]

png('out/diff.png',width=2000,height=2000)
par(mar=c(4,4,1,9))
plotDNA(apply(onlyDiff,1,paste,collapse='')[hiv$donor],groups=paste(hiv$select,hiv$Pair.ID)[hiv$donor])
dev.off()
png('out/align.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotDNA(align[hiv$seqId,'seq'],groups=paste(ifelse(hiv$donor,'D','R'),hiv$Pair.ID))
dev.off()


aa<-readFaDir('data','^[^LA].*aa.fasta.gz$') #don't read LTR since no AA or AlignSeq
bigAA<-tapply(aa$seq,aa$name,paste,collapse='')
protIndices<-withAs(xx=aa[aa$name==refName,],tapply(xx$seq,xx$file,function(x)cumsum(strsplit(x,'')[[1]]!='-')))
protIndices<-sprintf('%s.%d',rep(sub('\\..*$','',names(protIndices)),sapply(protIndices,length)),unlist(protIndices)) protIndices<-sprintf('%s.%s',protIndices,ave(protIndices,protIndices,FUN=function(x)1:length(x)))
aaMat<-seqSplit(bigAA)
colnames(aaMat)<-protIndices
rownames(aaMat)<-names(bigAA)
#refPos<-cumsum(aaMat[refName,]!='-')
aaMat<-aaMat[hiv$seqId,]
#colnames(aaMat)<-sprintf('aa%s.%s',refPos,ave(refPos,refPos,FUN=function(x)1:length(x)))
nBase<-apply(aaMat,2,function(x)sum(x!='-'))
seqMat<-seqMat[,nBase>10]
onlyDiffAA<-aaMat[,apply(aaMat,2,function(x)length(unique(x[x!='-']))>1)]


png('out/diffAA.png',width=2000,height=2000)
par(mar=c(4,4,1,9))
plotAA(apply(onlyDiffAA,1,paste,collapse='')[hiv$donor],groups=paste(hiv$select,hiv$Pair.ID)[hiv$donor])
dev.off()
png('out/alignAA.png',width=2000,height=2000)
par(mar=c(4,4,1,6))
plotAA(apply(aaMat,1,paste,collapse=''),groups=paste(ifelse(hiv$donor,'D','R'),hiv$Pair.ID))
dev.off()


png('out/diffAA_ifnb',width=2000,height=2000,res=200)
par(mar=c(4,4,1,9))
plotAA(apply(onlyDiffAA,1,paste,collapse='')[hiv$donor][order(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml[hiv$donor])])
dev.off()
png('out/diff_ifnb.png',width=2000,height=2000,res=200)
par(mar=c(4,4,1,9))
plotDNA(apply(onlyDiff,1,paste,collapse='')[hiv$donor][order(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml[hiv$donor])])
dev.off()
png('out/diffAA_ifna',width=2000,height=2000,res=200)
par(mar=c(4,4,1,9))
plotAA(apply(onlyDiffAA,1,paste,collapse='')[hiv$donor][order(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml[hiv$donor])])
dev.off()
png('out/diff_ifna.png',width=2000,height=2000,res=200)
par(mar=c(4,4,1,9))
plotDNA(apply(onlyDiff,1,paste,collapse='')[hiv$donor][order(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml[hiv$donor])])
dev.off()


