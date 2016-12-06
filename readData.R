library(levenR)
library(vipor)
library(dnar)
library(dnaplotr)
library(parallel)
library(xlsx)

#hiv<-read.csv('data/data.csv',stringsAsFactors=FALSE)
hiv<-read.xlsx("data/Final all data master 101116 copy 2.xlsx",sheetIndex=1,stringsAsFactors=FALSE)
hiv<-hiv[hiv$Renamed!='',]
hiv<-hiv[apply(hiv,1,function(x)!all(is.na(x))),apply(hiv,2,function(x)!all(is.na(x)))]
colnames(hiv)<-sub('\\.+$','',sub('^X\\.','',colnames(hiv)))
hiv<-hiv[hiv$Renamed!='CH742.SE.082708.BE.1',]
hiv$seq<-toupper(gsub('[ \n]','',hiv$Sequence))
hiv$select<-sapply(strsplit(hiv$Renamed,'\\.'),'[',4)
hiv$fluid<-sapply(strsplit(hiv$Renamed,'\\.'),'[',2)
hiv$donorRec<-sub(' ','-',hiv$Donor.or.Recipient)
hiv$sample<-sub('(-[0-9])(.*)$','\\2\\1',paste(hiv$donorRec,hiv$Pair.ID))
hiv$donor<-grepl('Donor',hiv$Donor.or.Recipient)
hiv$sampleSelect<-paste(hiv$donorRec,hiv$Pair.ID,hiv$select)
hiv$sampleFluid<-paste(hiv$donorRec,hiv$Pair.ID,hiv$fluid)
hiv$sampleFluidSelect<-paste(hiv$donorRec,hiv$Pair.ID,hiv$fluid,hiv$select)
hiv$fluidSelectDonor<-paste(ifelse(hiv$donor,'DO','RE'),ifelse(hiv$fluid=="PL",'PL','GE'),hiv$select)
hiv$fluidSelectDonor2<-paste(ifelse(hiv$donor,'Donor','Recipient'),ifelse(hiv$fluid=="PL",'Plasma','Genital'),ifelse(hiv$select=='UT','Untreated',ifelse(hiv$select=='A2','IFNA2','IFNB')))
hiv$seqId<-sprintf('seq%d',1:nrow(hiv))
hiv$baseName<-sub('\\..*$','',hiv$Renamed)
hiv$nameFluidSelect<-paste(hiv$baseName,hiv$fluid,hiv$select)
hiv$isGenital<-hiv$fluid!='PL'
donors<-with(hiv[hiv$donor,],tapply(baseName,Pair.ID,unique))
recs<-with(hiv[!hiv$donor,],tapply(baseName,Pair.ID,unique))
pairNames<-mapply(function(x,y)paste(paste(x,collapse='/'),paste(y,collapse='/'),sep='-'),donors,recs)
write.fa(hiv$seqId,hiv$seq,'hiv.fa')

#deal with column name case inconsistency
colnames(hiv)<-sub('\\.donor\\.','.Donor.',colnames(hiv))
hiv$meanIfna<-(hiv$IFNa2.Single.Donor.cells.IC50..pg..ml+hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml)/2

#calculate scaled replicative capacity
hiv[hiv$select=='UT','maxSD']<-ave(hiv[hiv$select=='UT','Replicative.capacity.Single.Donor.cells.p24.d7'], hiv[hiv$select=='UT','Pair.ID'], FUN=max)
hiv[hiv$select=='UT','maxPD']<-ave(hiv[hiv$select=='UT','Replicative.capacity.Pooled.Donor.cells.p24.d7'], hiv[hiv$select=='UT','Pair.ID'], FUN=max)
#fill in the treated ones
hiv$maxSD<-ave(hiv$maxSD,hiv$Pair.ID,FUN=function(x)max(x,na.rm=TRUE))
hiv$maxPD<-ave(hiv$maxPD,hiv$Pair.ID,FUN=function(x)max(x,na.rm=TRUE))
hiv$propSD<-hiv$Replicative.capacity.Single.Donor.cells.p24.d7/hiv$maxSD
hiv$propPD<-hiv$Replicative.capacity.Pooled.Donor.cells.p24.d7/hiv$maxPD
#take mean of pooled and single
hiv$meanRepCap<-(hiv$propSD+hiv$propPD)/2
write.csv(hiv[,c('Renamed','Original.name','Pair.ID','select','Replicative.capacity.Single.Donor.cells.p24.d7','Replicative.capacity.Pooled.Donor.cells.p24.d7','maxSD','maxPD','propSD','propPD','meanRepCap')],'out/repCap.csv')

source('readVarGlyco.R')

#Read in precise vres values to figure out which are censored
vresPrecise<-read.csv('data/Vres data for scott v2.csv',stringsAsFactors=FALSE)
vresPrecise<-vresPrecise[apply(is.na(vresPrecise),1,sum)==0,]
rownames(vresPrecise)<-vresPrecise$Renamed
vresPrecise$isCensor<-vresPrecise$IFN.beta.Pooled.Donor.p24.at.0.44.pg.ml==.1
vresPrecise$vres<-vresPrecise$IFN.beta.Pooled.Donor.p24.at.0.44.pg.ml/vresPrecise$IFN.beta.expt.Replicative.capacity..p24...d7.ng.ml*100
hiv$vres<-vresPrecise[hiv$Renamed,'vres']
hiv$vresCensor<-vresPrecise[hiv$Renamed,'isCensor']
if(any(is.na(hiv$vres)))stop(simpleError('Unknown vres censoring'))
if(any(is.na(hiv$vresCensor)))stop(simpleError('Unknown vres censoring'))

hiv$minVres<-.1/hiv$IFN.beta.expt.Replicative.capacity..p24...d7.ng.ml*100
hiv$isMinVres<-round(hiv$vres,2)<=round(hiv$minVres,2)
if(any(round(hiv$vres,5)<round(hiv$minVres)))stop(simpleError('Vres < min Vres'))
if(any(round(hiv$vres,5)<=round(hiv$minVres,5)&!hiv$vresCensor))stop(simpleError('Vres == min Vres and not censored'))
if(any(round(hiv$vres,5)>round(hiv$minVres,5)*1.02&hiv$vresCensor))stop(simpleError('Vres > min Vres and censored'))


targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Replicative capacity (p24)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml'='IFNa2 IC50 (pg/ml)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml'='IFNbeta IC50 (pg/ml)',
  'Residual.Pooled.Donor.cells..1500U..UT'='IFNa2 Vres',
  'vres'='IFNbeta Vres',
  'p24.release.No.IFN'='p24 release'
  #'p24.release.With.IFNa..500.U.ml'='p24 release with IFNa2',
  #'Autologous.IC50'='Autologous IC50',
  #'Bnaber.IC50'='Bnaber IC50'
)
targetColTransform<-structure(rep('log',length(targetCols)),names=names(targetCols))
#targetColTransform['Replicative.capacity.Pooled.Donor.cells.p24.d7']<-'identity'
targetColTransform[grep('p24.release',names(targetColTransform))]<-'logit'
targetColTransform[names(targetColTransform)=='vres']<-'logit'
targetColTransform[names(targetColTransform)=='Residual.Pooled.Donor.cells..1500U..UT']<-'logit'
targetColCensorDown<-rep(list(c()),length(targetCols))
names(targetColCensorDown)<-names(targetCols)
targetColCensorDown[['Autologous.IC50']]<-rep(20,nrow(hiv))
targetColCensorDown[['Bnaber.IC50']]<-rep(20,nrow(hiv))
targetColCensorDown[['vres']]<-hiv$minVres
goodTargetCols<-targetCols[1:7]

#note using log10 instead of log to make plotting easier
logit<-function(p)log10(p)-log10(1-p)
invLogit<-function(x)10^(x)/(10^(x)+1)

pairColors<-c('3'='#999999','4'='#99CC33','7'='#CC6699','1'='#9999CC','2'='#99CCCC','5'='#CC9966','6'='#FF9966')

#output final csv
desiredCols<-c('Name'='Renamed','Subtype'='Subtype','Gender'='Gender','Donor/Recipient'='Donor.or.Recipient','Pair ID'='Pair.ID','Fluid'='fluid','Selection'='select',structure(names(targetCols),.Names=targetCols),'Censored IFNbeta Vres'='isMinVres')
finalCsv<-hiv[,desiredCols]
names(finalCsv)<-names(desiredCols)
finalCsv$Name[finalCsv$Fluid!='PL']<-sub('\\.UT\\.','.',finalCsv$Name[finalCsv$Fluid!='PL'])
write.csv(finalCsv,'out/Iyer2016_Data.csv',row.names=FALSE)



