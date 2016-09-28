library('rstan')
library('dnar')
library('vioplot')

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
    int<lower=1> groupIds[N];
    int<lower=0> nPair;
    int<lower=1> pairIds[N];
    int<lower=0> nRecipient;
    int<lower=0> recipientIds[N];
    int<lower=0> nGroupTypes;
    int<lower=1> groupTypes[nGroup];
    int<lower=0> nGenital;
    int<lower=0> nCladeB;
    int<lower=1> cladeBId[N];
  }
  parameters {
    real metaDonorMu;
    real<lower=0> metaDonorSd;
    real metaRecipientMu;
    real<lower=0> metaRecipientSd;
    real metaGenitalMu;
    real<lower=0> metaGenitalSd;
    real metaCladeMu;
    real<lower=0> metaCladeSd;
    real donors[nPair];
    real<lower=0> sigmas[nGroup];
    real<lower=0> metaSigmaMu[nGroupTypes];
    real<lower=0> metaSigmaSigma[nGroupTypes];
    real genitals[nGenital];
    real recipients[nRecipient];
    real clades[nCladeB];
  }
  transformed parameters{
    real indivMu[N];
    real<lower=0> indivSigma[N];
    for (ii in 1:N){
      indivMu[ii] = metaDonorMu+donors[pairIds[ii]]*metaDonorSd;
      if(isRecipient[ii])indivMu[ii]=indivMu[ii]+metaRecipientMu+recipients[recipientIds[ii]]*metaRecipientSd;
      if(isGenital[ii])indivMu[ii]=indivMu[ii]+metaGenitalMu+genitals[pairIds[ii]]*metaGenitalSd;
      if(isCladeB[ii] && isRecipient[ii])indivMu[ii]=indivMu[ii]+metaCladeMu+clades[cladeBId[ii]]*metaCladeSd;
      indivSigma[ii]=metaSigmaMu[groupTypes[groupIds[ii]]]+sigmas[groupTypes[groupIds[ii]]]*metaSigmaSigma[groupTypes[groupIds[ii]]];
    }
  }
  model {
    metaSigmaMu ~ gamma(2,.5);
    metaSigmaSigma ~ gamma(2,.5);
    metaDonorSd ~ gamma(2,.5);
    metaRecipientSd ~ gamma(2,.5);
    metaGenitalSd ~ gamma(2,.5);
    metaCladeSd ~ gamma(2,.5);
    donors ~ normal(0,1);
    recipients ~ normal(0,1);
    genitals ~ normal(0,1);
    clades ~ normal(0,1);
    sigmas ~ normal(0,1);
    ic50 ~ normal(indivMu,indivSigma);
  }
"

targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Pooled donor\nReplicative capacity (day 7 p24)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml.'='IFNbeta IC50 (pg/ml)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'='IFNa2 IC50 (pg/ml)'
)
fits<-lapply(names(targetCols),function(targetCol){
  groupTypes<-sapply(1:max(hiv$group),function(zz)paste(ifelse(hiv[hiv$group==zz,'fluid'][1]=='PL','PL','GE'),ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec')))
  #just group by donor or recipient
  #groupTypes<-sapply(1:max(hiv$group),function(zz)ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec'))
  cladeBs<-unique(hiv$Pair.ID..[hiv$Subtype=='B'])
  notCladeBs<-unique(hiv$Pair.ID..[hiv$Subtype!='B'])
  #note 99999 is a arbitrarily high number for non clade Bs (should never be called within Stan due to if(cladeB))
  cladeBIds<-structure(c(1:length(cladeBs),rep(99999,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  recipients<-unique(hiv[!hiv$donor,'sample'][order(hiv$Pair.ID..[!hiv$donor])])
  notRecipients<-unique(hiv[hiv$donor,'sample'])
  recipientIds<-structure(c(1:length(recipients),rep(99999,length(notRecipients))),.Names=c(recipients,notRecipients))
  dat<-withAs('xx'=hiv[hiv$select=='UT'&!is.na(hiv[,targetCol]),],list(
    ic50=log10(xx[,targetCol]),
    N=nrow(xx),
    isCladeB=as.integer(xx$Subtype=='B'),
    isGenital=as.integer(xx$fluid!='PL'),
    nGenital=max(xx[xx$fluid!='PL','Pair.ID..']),
    isRecipient=as.integer(!xx$donor),
    nGroup=max(xx$group),
    groupIds=xx$group,
    nPair=max(xx$Pair.ID..),
    pairIds=xx$Pair.ID..,
    nRecipient=length(unique(xx[!xx$donor,'sample'])),
    recipientIds=recipientIds[xx$sample],
    nGroupTypes=length(unique(groupTypes)),
    groupTypes=as.numeric(as.factor(groupTypes)),
    nCladeB=length(cladeBs),
    cladeBId=cladeBIds[as.character(xx$Pair.ID..)]
  ))
  fit <- cacheOperation(sprintf('work/stan%s.Rdat',targetCol),stan,model_code = stanCode, data = dat, iter=50000, chains=nThreads,thin=25,control=list(adapt_delta=.99))
  return(list('fit'=fit,'dat'=dat))
})
names(fits)<-names(targetCols)

for(targetCol in names(targetCols)){
  fit<-fits[[targetCol]][['fit']]
  dat<-fits[[targetCol]][['dat']]
  #
  sims<-as.array(fit)
  dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
  colnames(sims)<-dimnames(as.array(fit))[[3]]
  #
  allPars<-c("metaDonorMu", "metaDonorSd", "metaRecipientMu", "metaRecipientSd", "metaGenitalMu", "metaGenitalSd","metaCladeMu","metaCladeSd","donors", "sigmas", "metaSigmaMu", "metaSigmaSigma", "genitals", "recipients", "clades")
  indivMuCols<-sprintf('indivMu[%d]',1:dat[['N']])
  indivSdCols<-sprintf('indivSigma[%d]',1:dat[['N']])
  simFits<-mapply(function(mu,sigma)rnorm(nrow(sims),sims[,mu],sims[,sigma]),indivMuCols,indivSdCols,SIMPLIFY=FALSE)
  #
  pdf(sprintf('out/bayes/bayesFit%s.pdf',targetCol),width=20,height=20)
    print(plot(fit,pars=allPars))
    print(traceplot(fit,pars=allPars))
    par(mar=c(4,5,2,0))
    plot(1,1,type='n',xlim=c(1,dat[['N']])+c(-1,1),ylim=range(unlist(simFits)),xlab='Virus ID',ylab=sprintf('Log10 %s',targetCols[targetCol]),xaxs='i',cex.axis=1.5,cex.lab=2,main='Posterior predictive distributions',cex.main=2)
    do.call(vioplot,c('x'=list(simFits[[1]]),simFits[-1],'add'=list(TRUE),'colMed'=list(NA),'col'=list('#00000033')))
    #points(rep(1:length(simFits),sapply(simFits,length)),unlist(simFits),pch='.',col='#00000033')
    points(1:length(dat[['ic50']]),dat[['ic50']],col='red',pch='+',lwd=2)
  dev.off()
  #
  indivRecipBeta<-sims[,grep('recipients\\[[0-9]\\]',colnames(sims))]
  metaBeta<-sims[,'metaRecipientMu']
  indivGenital<-sims[,grep('genitals\\[[0-9]\\]',colnames(sims))]
  metaGenital<-sims[,'metaGenitalMu']
  indivClade<-sims[,grep('clades\\[[0-9]\\]',colnames(sims))]
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
  pdf(sprintf('out/bayes/bayes%s.pdf',targetCol))
    par(mfrow=c(2,1),las=1,mar=c(3,3.6,1.1,.1))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.7,0),log='x',main='Individuals',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    meanBin<-(bins[-length(bins)]+bins[-1])/2
    apply(indivTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(indivGenitalTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(indivCladeTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    abline(v=1,lty=2)
    logAxis(1,mgp=c(3,.5,0))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Population',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaTabs,0,0),col='#0000FF44',border='#0000FF99')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,genitalTabs,0,0),col='#FF000044',border='#FF000099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,cladeTabs,0,0),col='#00FF0044',border='#00FF0099')
    abline(v=1,lty=2)
    legend('topleft',c('Recipient','Clade B','Genital'),fill=c('#0000FF44','#00FF0044','#FF000044'),border=c('#0000FF99','#00FF0099','#FF000099'),inset=.02)
    logAxis(1,mgp=c(3,.5,0))
  dev.off()
}


