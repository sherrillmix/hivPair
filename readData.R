library(levenR)
library(vipor)
library(dnar)
library(dnaplotr)
library(parallel)
#library(beeswarm)

hiv<-read.csv('data/data.csv',stringsAsFactors=FALSE)
hiv<-hiv[hiv$Renamed!='',]
hiv<-hiv[,apply(hiv,2,function(x)!all(is.na(x)))]
hiv$seq<-toupper(gsub('[ \n]','',hiv$Sequence))
hiv$select<-sapply(strsplit(hiv$Renamed,'\\.'),'[',4)
hiv$fluid<-sapply(strsplit(hiv$Renamed,'\\.'),'[',2)
hiv$donorRec<-sub(' ','-',hiv$Donor.or.Recipient)
hiv$sample<-paste(hiv$donorRec,hiv$Pair.ID..)
hiv$donor<-grepl('Donor',hiv$Donor.or.Recipient)
hiv$sampleSelect<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$select)
hiv$sampleFluid<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$fluid)
hiv$sampleFluidSelect<-paste(hiv$donorRec,hiv$Pair.ID..,hiv$fluid,hiv$select)
hiv$fluidSelectDonor<-paste(ifelse(hiv$donor,'DO','RE'),ifelse(hiv$fluid=="PL",'PL','GE'),hiv$select)
hiv$seqId<-sprintf('seq%d',1:nrow(hiv))
hiv$baseName<-sub('\\..*$','',hiv$Renamed)
hiv$isGenital<-hiv$fluid!='PL'
donors<-with(hiv[hiv$donor,],tapply(baseName,Pair.ID..,unique))
recs<-with(hiv[!hiv$donor,],tapply(baseName,Pair.ID..,unique))
pairNames<-mapply(function(x,y)paste(paste(x,collapse='/'),paste(y,collapse='/'),sep='-'),donors,recs)
write.fa(hiv$seqId,hiv$seq,'hiv.fa')

#deal with column name case inconsistency
colnames(hiv)<-sub('\\.donor\\.','.Donor.',colnames(hiv))
hiv$meanIfna<-(hiv$IFNa2.Single.Donor.cells.IC50..pg..ml.+hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.)/2

#calculate scaled replicative capacity
hiv[hiv$select=='UT','maxSD']<-ave(hiv[hiv$select=='UT','Replicative.capacity.Single.Donor.cells.p24.d7'], hiv[hiv$select=='UT','Pair.ID..'], FUN=max)
hiv[hiv$select=='UT','maxPD']<-ave(hiv[hiv$select=='UT','Replicative.capacity.Pooled.Donor.cells.p24.d7'], hiv[hiv$select=='UT','Pair.ID..'], FUN=max)
#fill in the treated ones
hiv$maxSD<-ave(hiv$maxSD,hiv$Pair.ID..,FUN=function(x)max(x,na.rm=TRUE))
hiv$maxPD<-ave(hiv$maxPD,hiv$Pair.ID..,FUN=function(x)max(x,na.rm=TRUE))
hiv$propSD<-hiv$Replicative.capacity.Single.Donor.cells.p24.d7/hiv$maxSD
hiv$propPD<-hiv$Replicative.capacity.Pooled.Donor.cells.p24.d7/hiv$maxPD
#take mean of pooled and single
hiv$meanRepCap<-(hiv$propSD+hiv$propPD)/2
write.csv(hiv[,c('Renamed','Original.name','Pair.ID..','select','Replicative.capacity.Single.Donor.cells.p24.d7','Replicative.capacity.Pooled.Donor.cells.p24.d7','maxSD','maxPD','propSD','propPD','meanRepCap')],'out/repCap.csv')


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

source('readVarGlyco.R')

targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Pooled donor\nReplicative capacity (day 7 p24)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml.'='IFNbeta IC50 (pg/ml)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'='IFNa2 IC50 (pg/ml)',
  'p24.release.With.IFNa..500.U.ml....'='p24 release with IFNa2',
  'p24.release.No.IFN....'='p24 release without IFNa2',
  'Autologous.IC50'='Autologous IC50',
  'Bnaber.IC50'='Bnaber IC50'
)
targetColTransform<-structure(rep('log',length(targetCols)),names=names(targetCols))
targetColTransform['Replicative.capacity.Pooled.Donor.cells.p24.d7']<-'identity'
targetColTransform[grep('p24.release.',names(targetColTransform))]<-'logit'
targetColCensorDown<-structure(rep(NA,length(targetCols)),names=names(targetCols))
targetColCensorDown['Autologous.IC50']<-20
targetColCensorDown['Bnaber.IC50']<-20
goodTargetCols<-targetCols[1:5]


