library('rstan')
library('dnar')

#for parallel
nThreads<-50
if(parallel::detectCores()>10){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

if(!exists('hiv'))source('readData.R')

stanCode<-"
  data {
    int<lower=0> N;
    real ic50[N];
    int<lower=0,upper=1> isCladeB[N];
    int<lower=0,upper=1> isGenital[N];
    int<lower=0,upper=1> isRecipient[N];
    int<lower=0> nGroup;
    int<lower=1> group[N];
    int<lower=0> nPair;
    int<lower=1> pair[N];
    int<lower=0> nGroupTypes;
    int<lower=1> groupTypes[nGroup];
    int<lower=0> nGenital;
    int<lower=0> nCladeB;
    int<lower=1> cladeBId[N];
  }
  parameters {
    real metaDonorMu;
    real<lower=0.0000001> metaDonorSd;
    real metaRecipientMu;
    real<lower=0.0000001> metaRecipientSd;
    real metaGenitalMu;
    real<lower=0.0000001> metaGenitalSd;
    real metaCladeMu;
    real<lower=0.0000001> metaCladeSd;
    real donors[nPair];
    real<lower=0.0000001> sigmas[nGroup];
    real<lower=0.0000001> metaSigmaSd[nGroupTypes];
    real metaSigmaMu[nGroupTypes];
    real genitals[nGenital];
    real recipients[nPair];
    real clade[nCladeB];
  }
  transformed parameters{
    real indivMu[N];
    for (ii in 1:N){
      indivMu[ii] = donors[pair[ii]]+recipients[pair[ii]]*isRecipient[ii];
      if(isGenital[ii])indivMu[ii]=indivMu[ii]+genitals[pair[ii]];
      if(isCladeB[ii])indivMu[ii]=indivMu[ii]+isRecipient[ii]*clade[cladeBId[ii]];
    }
  }
  model {
    donors ~ normal(metaDonorMu,metaDonorSd);
    recipients ~ normal(metaRecipientMu,metaRecipientSd);
    genitals ~ normal(metaGenitalMu,metaGenitalSd);
    for(ii in 1:nGroup)sigmas[ii] ~ normal(metaSigmaMu[groupTypes[ii]],metaSigmaSd[groupTypes[ii]]);
    for (ii in 1:N)ic50[ii] ~ normal(indivMu[ii],sigmas[group[ii]]);
  }
"

targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Pooled donor\nReplicative capacity (day 7 p24)',
  'meanRepCap'='Mean replicative capacity\n(proportion of maximum day 7 p24)',
  'meanIfna'='IFNa2 IC50 (pg/ml)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml.'='IFNbeta IC50 (pg/ml)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'='IFNa2 IC50 (pg/ml)'
)
fits<-lapply(targetCols,function(targetCol){
  groupTypes<-sapply(1:max(hiv$group),function(zz)paste(ifelse(hiv[hiv$group==zz,'fluid'][1]=='PL','PL','GE'),ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec')))
  cladeBs<-unique(hiv$Pair.ID..[hiv$Subtype=='B'])
  notCladeBs<-unique(hiv$Pair.ID..[hiv$Subtype!='B'])
  #note 99999 is a arbitrarily high number for non clade Bs (should never be called within Stan due to if(cladeB))
  cladeBIds<-structure(c(1:length(cladeBs),rep(99999,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  dat<-withAs('xx'=hiv[hiv$select=='UT'&!is.na(hiv[,targetCol]),],list(
    ic50=log10(xx[,targetCol]),
    N=nrow(xx),
    isCladeB=as.integer(xx$Subtype=='B'),
    isGenital=as.integer(xx$fluid!='PL'),
    nGenital=max(xx[xx$fluid!='PL','Pair.ID..']),
    isRecipient=as.integer(!xx$donor),
    nGroup=max(xx$group),
    group=xx$group,
    pair=xx$Pair.ID..,
    nPair=max(xx$Pair.ID..),
    nGroupTypes=length(unique(groupTypes)),
    groupTypes=as.numeric(as.factor(groupTypes)),
    nCladeB=length(cladeBs),
    cladeBId=cladeBIds[as.character(xx$Pair.ID..)]
  ))
  fit <- cacheOperation(sprintf('work/stan%s.Rdat',targetCol),stan,model_code = stanCode, data = dat, iter=30000, chains=nThreads,thin=20)
  return(fit)
})
names(fits)<-targetCols

for(targetCol in targetCols){
  fit<-fits[[targetCol]]
  allPars<-c("metaDonorMu", "metaDonorSd", "metaRecipientMu", "metaRecipientSd", "metaGenitalMu", "metaGenitalSd","metaCladeMu","metaCladeSd","donors", "sigmas", "metaSigmaSd", "metaSigmaMu", "genitals", "recipients", "clade")
  pdf(sprintf('out/bayesFit%s.pdf',targetCol),width=20,height=20)
    print(plot(fit,pars=allPars))
    print(traceplot(fit,pars=allPars))
  dev.off()
  #
  sims<-as.array(fit)
  dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
  colnames(sims)<-dimnames(as.array(fit))[[3]]
  #
  indivRecipBeta<-sims[,grep('recipients\\[[0-9]\\]',colnames(sims))]
  metaBeta<-sims[,'metaRecipientMu']
  indivGenital<-sims[,grep('genitals\\[[0-9]\\]',colnames(sims))]
  metaGenital<-sims[,'metaGenitalMu']
  indivClade<-sims[,grep('clade\\[[0-9]\\]',colnames(sims))]
  metaClade<-sims[,'metaCladeMu']
  xlim<-c(-.7,2.3)
  #
  bins<-seq(xlim[1],xlim[2],.01)
  indivTabs<-apply(indivRecipBeta,2,function(x)table(cut(x,bins))/length(x))
  metaTabs<-table(cut(metaBeta,bins))/length(metaBeta)
  indivGenitalTabs<-apply(indivGenital,2,function(x)table(cut(x,bins))/length(x))
  genitalTabs<-table(cut(metaGenital,bins))/length(metaGenital)
  indivCladeTab<-apply(indivClade,2,function(x)table(cut(x,bins))/length(x))
  cladeTabs<-table(cut(metaClade,bins))/length(metaClade)
  #
  pdf(sprintf('out/bayes%s.pdf',targetCol))
    par(mfrow=c(2,1),las=1,mar=c(4,3.6,1.1,.1))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log='x',main='Individuals',xaxs='i',main=names(targetCols)[targetCols==targetCol])
    #axis(1,prettyLabels,sapply(prettyLabels,function(x)as.expression(bquote(2^.(x)))),las=1)
    #axis(1,log2(10^prettyLabels),ifelse(prettyLabels==0,1,sapply(prettyLabels,function(x)as.expression(bquote(10^.(x))))),las=1)
    title(xlab='Fold increase',mgp=c(1.5,1,0))
    meanBin<-(bins[-length(bins)]+bins[-1])/2
    apply(indivTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(indivGenitalTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(indivCladeTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    abline(v=1,lty=2)
    logAxis(1)
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Population',xaxs='i')
    title(xlab='Fold increase',mgp=c(2,1,0))
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaTabs,0,0),col='#0000FF44',border='#0000FF99')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,genitalTabs,0,0),col='#FF000044',border='#FF000099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,cladeTabs,0,0),col='#00FF0044',border='#00FF0099')
    abline(v=1,lty=2)
    legend('topleft',c('Recipient','Clade B','Genital'),fill=c('#0000FF44','#00FF0044','#FF000044'),border=c('#0000FF99','#00FF0099','#FF000099'),inset=.02)
    logAxis(1)
  dev.off()
}


