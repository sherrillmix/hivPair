library('rstan')

#for parallel
if(parallel::detectCores()>10){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  nThreads=20
}else{
  nThreads=4
}

if(!exists('hiv'))source('readData.R')

targetCol<-'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'
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
  }
  parameters {
    real metaMuMu;
    real<lower=0.0000001> metaMuSd;
    real metaRecipientMu;
    real<lower=0.0000001> metaRecipientSd;
    real metaGenitalMu;
    real<lower=0.0000001> metaGenitalSd;
    real mus[nPair];
    real<lower=0.0000001> sigmas[nGroup];
    real<lower=0.0000001> metaSigmaSd[nGroupTypes];
    real metaSigmaMu[nGroupTypes];
    real genitals[nGenital];
    real recipient[nPair];
    real clade;
  }
  transformed parameters{
    real indivMu[N];
    for (ii in 1:N){
      indivMu[ii] = mus[pair[ii]]+recipient[pair[ii]]*isRecipient[ii]+clade*isCladeB[ii];
      if(isGenital[ii])indivMu[ii]=indivMu[ii]+genitals[pair[ii]];
    }
  }
  model {
    mus ~ normal(metaMuMu,metaMuSd);
    recipient ~ normal(metaRecipientMu,metaRecipientSd);
    genitals ~ normal(metaGenitalMu,metaGenitalSd);
    for(ii in 1:nGroup)sigmas[ii]~normal(metaSigmaMu[groupTypes[ii]],metaSigmaSd[groupTypes[ii]]);
    for (ii in 1:N){
      ic50[ii]~normal(
        indivMu[ii],
        sigmas[group[ii]]
      );
    }
  }
"

groupTypes<-sapply(1:max(hiv$group),function(zz)paste(ifelse(hiv[hiv$group==zz,'fluid'][1]=='PL','PL','GE'),ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec')))
dat<-withAs('xx'=hiv[hiv$select=='UT',],list(
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
  groupTypes=as.numeric(as.factor(groupTypes))
))
fit <- stan(model_code = stanCode, data = dat, iter=10000, chains=nThreads)

allPars<-c("metaMuMu", "metaMuSd", "metaRecipientMu", "metaRecipientSd", "metaGenitalMu", "metaGenitalSd", "mus", "sigmas", "metaSigmaSd", "metaSigmaMu", "genitals", "recipient", "clade")
pdf('test.pdf',width=20,height=20)
  print(plot(fit,pars=allPars))
  print(traceplot(fit,pars=allPars))
dev.off()

sims<-as.array(fit)
dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
colnames(sims)<-dimnames(as.array(fit))[[3]]

indivRecipBeta<-sims[,grep('recipient\\[[0-9]\\]',colnames(sims))]
metaBeta<-sims[,'metaRecipientMu']


xlim<-c(-.1,1.2)
#prettyLabels<-pretty(c(-5,25))
bins<-seq(xlim[1],xlim[2],.01)
indivTabs<-apply(indivRecipBeta,2,function(x)table(cut(x,bins))/length(x))
metaTabs<-table(cut(metaBeta,bins))/length(metaBeta)

pdf('test.pdf')
  par(las=1,mar=c(4,3.6,1.1,.1))
  plot(1,1,type='n',xlim=xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0))
  #axis(1,prettyLabels,sapply(prettyLabels,function(x)as.expression(bquote(2^.(x)))),las=1)
  #axis(1,log2(10^prettyLabels),ifelse(prettyLabels==0,1,sapply(prettyLabels,function(x)as.expression(bquote(10^.(x))))),las=1)
  #title(xlab='Fold increase in DNA over air swabs',mgp=c(2,1,0))
  meanBin<-(bins[-length(bins)]+bins[-1])/2
  apply(indivTabs,2,function(xx)polygon(c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11'))
  polygon(c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaTabs,0,0),col='#FF000044')
  abline(v=0,lty=2)
dev.off()


