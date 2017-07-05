if(!exists('hiv'))source("readData.R")
refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
align<-read.fa('data/AlignSeq.nt.fasta.gz')
rownames(align)<-align$name
seqMat<-seqSplit(align$seq)
rownames(seqMat)<-align$name
refPos<-cumsum(seqMat[refName,]!='-')
seqMat<-seqMat[hiv$seqId,]
if(any(degap(apply(seqMat,1,paste,collapse=''))!=hiv$seq))stop(simpleError('hiv metadata and alignment mismatched'))
colnames(seqMat)<-sprintf('nt%s.%s',refPos,ave(refPos,refPos,FUN=function(x)1:length(x)))
nBase<-apply(seqMat,2,function(x)sum(x!='-')) startEnd<-range(which(nBase/nrow(seqMat)>.9))
seqMat<-seqMat[,nBase>10&1:ncol(seqMat) %in% startEnd[1]:startEnd[2]]
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
protIndices<-sprintf('%s.%d',rep(sub('\\..*$','',names(protIndices)),sapply(protIndices,length)),unlist(protIndices))
protIndices<-sprintf('%s.%s',protIndices,ave(protIndices,protIndices,FUN=function(x)1:length(x)))
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


allNonGap<-range(which(apply(seqMat!='-',2,all)))
noGapSeqs<-apply(seqMat[,allNonGap[1]:allNonGap[2]],1,paste,collapse='')
xx<-degap(noGapSeqs)
gcPos<-gregexpr('CG',xx)
gcs<-sapply(gcPos,length)/nchar(xx)
gcaPos<-gregexpr('[AT]CG[AT]',xx)
gcas<-sapply(gcaPos,length)/nchar(xx)
gcaPerGc<-gcas/gcs
pdf('out/gc.pdf')
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcs[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='CG/nucleotide',las=1)
  for(ii in unique(hiv$Pair.ID)){
    main<-paste(unique(hiv$baseName[hiv$Pair.ID==ii][order(!hiv$donor[hiv$Pair.ID==ii])]),collapse=',')
    vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$Pair.ID==ii&hiv$select=='UT'],gcs[hiv$Pair.ID==ii&hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID][hiv$Pair.ID==ii],ylab='CG/nucleotide',main=main,las=1,ylim=range(gcs))
  }
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcas[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='[AT]CG[AT]/nucleotide',las=1)
  for(ii in unique(hiv$Pair.ID)){
    main<-paste(unique(hiv$baseName[hiv$Pair.ID==ii][order(!hiv$donor[hiv$Pair.ID==ii])]),collapse=',')
    vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$Pair.ID==ii&hiv$select=='UT'],gcas[hiv$Pair.ID==ii&hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID][hiv$Pair.ID==ii],ylab='[AT]CG[AT]/nucleotide',main=main,las=1,ylim=range(gcas))
  }
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcaPerGc[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='[AT]CG[AT]/CG',las=1)
  plotDNA(noGapSeqs)
dev.off()
