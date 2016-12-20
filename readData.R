#selection abbreviations
selectionExpand<-c('UT'='Untreated','A2'='IFNA2','BE'='IFNB')

#read data
hiv<-read.csv('data/Iyer2016_Data.csv',stringsAsFactors=FALSE,check.names=FALSE)
hiv$baseName<-sub('\\..*$','',hiv$Name)
hiv$isGenital<-hiv$Fluid!='PL'
hiv$isDonor<-hiv[,'Donor/Recipient']=='Donor'
hiv$fluidSelectDonor<-paste(ifelse(hiv$isDonor,'Donor','Recipient'),ifelse(hiv$isGenital,'Genital','Plasma'),selectionExpand[hiv$Selection])
hiv$sample<-paste(hiv[,'Donor/Recipient'],hiv[,'Pair ID'])
hiv$sampleFluidSelect<-paste(hiv$sample,hiv[,'Fluid'],hiv[,'Selection'])
hiv$sampleSelect<-paste(sub(' ','-',hiv[,'Donor/Recipient']),hiv[,'Pair ID'],hiv$Selection)

#columns with data
targetCols<-c(
  "Env/RT ratio",
  "Infectivity (RLU/pg RT)",
  "Replicative capacity (ng p24/ml)",
  "IFNa2 IC50 (pg/ml)",
  "IFNbeta IC50 (pg/ml)",
  "IFNa2 Vres",
  "IFNbeta Vres",
  "p24 release"
)
goodTargetCols<-targetCols[apply(is.na(hiv[,targetCols]),2,mean)<.1]

#transformations for each variable
targetColTransform<-structure(rep('log',length(targetCols)),names=targetCols)
targetColTransform[grep('p24 release|Vres',names(targetColTransform))]<-'logit'

#censored observations
targetColCensorDown<-rep(list(c()),length(targetCols))
names(targetColCensorDown)<-targetCols
targetColCensorDown[['IFNbeta Vres']]<-ifelse(hiv[,'Censored IFNbeta Vres'],hiv[,'IFNbeta Vres'],0)

if(!dir.exists('out'))dir.create('out')

pairColors<-c('1'='#999999','2'='#99CC33','3'='#CC6699','4'='#9999CC','5'='#99CCCC','6'='#CC9966','7'='#FF9966')

#note using log10 instead of log to make plotting easier
logit<-function(p)log10(p)-log10(1-p)
invLogit<-function(x)10^(x)/(10^(x)+1)
