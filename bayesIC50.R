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
    int<lower=0> nAlpha;
    int<lower=1> alphaIds[N];
    int<lower=0> nBeta;
    int<lower=1> betaIds[N];
    int<lower=0> nRecipientAlpha;
    int<lower=1> recipientAlphaIds[N];
    int<lower=0> nRecipientBeta;
    int<lower=1> recipientBetaIds[N];
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
    real metaAlphaMu;
    real<lower=0> metaAlphaSd;
    real metaBetaMu;
    real<lower=0> metaBetaSd;
    real metaRecipientAlphaMu;
    real<lower=0> metaRecipientAlphaSd;
    real metaRecipientBetaMu;
    real<lower=0> metaRecipientBetaSd;
    real alphas[nAlpha];
    real betas[nBeta];
    real recipientAlphas[nRecipientAlpha];
    real recipientBetas[nRecipientBeta];
  }
  transformed parameters{
    real indivMu[N];
    real<lower=0> indivSigma[N];
    for (ii in 1:N){
      indivMu[ii] = metaDonorMu+donors[pairIds[ii]]*metaDonorSd;
      if(isRecipient[ii])indivMu[ii]=indivMu[ii]+metaRecipientMu+recipients[recipientIds[ii]]*metaRecipientSd;
      if(isGenital[ii])indivMu[ii]=indivMu[ii]+metaGenitalMu+genitals[pairIds[ii]]*metaGenitalSd;
      if(isCladeB[ii] && isRecipient[ii])indivMu[ii]=indivMu[ii]+metaCladeMu+clades[cladeBId[ii]]*metaCladeSd;
      if(recipientAlphaIds[ii]<999)indivMu[ii]=indivMu[ii]+metaRecipientAlphaMu+recipientAlphas[recipientAlphaIds[ii]]*metaRecipientAlphaSd;
      if(recipientBetaIds[ii]<999)indivMu[ii]=indivMu[ii]+metaRecipientBetaMu+recipientBetas[recipientBetaIds[ii]]*metaRecipientBetaSd;
      if(alphaIds[ii]<999)indivMu[ii]=indivMu[ii]+metaAlphaMu+alphas[alphaIds[ii]]*metaAlphaSd;
      if(betaIds[ii]<999)indivMu[ii]=indivMu[ii]+metaBetaMu+betas[betaIds[ii]]*metaBetaSd;
      indivSigma[ii]=metaSigmaMu[groupTypes[groupIds[ii]]]+sigmas[groupIds[ii]]*metaSigmaSigma[groupTypes[groupIds[ii]]];
    }
  }
  model {
    metaSigmaMu ~ gamma(1,2);
    metaSigmaSigma ~ gamma(1,2);
    metaDonorSd ~ gamma(1,2);
    metaRecipientSd ~ gamma(1,2);
    metaGenitalSd ~ gamma(1,2);
    metaCladeSd ~ gamma(1,2);
    metaAlphaSd ~ gamma(1,2);
    metaBetaSd ~ gamma(1,2);
    metaRecipientAlphaSd ~ gamma(1,2);
    metaRecipientBetaSd ~ gamma(1,2);
    donors ~ normal(0,1);
    recipients ~ normal(0,1);
    genitals ~ normal(0,1);
    clades ~ normal(0,1);
    sigmas ~ normal(0,1);
    alphas ~ normal(0,1);
    betas ~ normal(0,1);
    recipientAlphas ~ normal(0,1);
    recipientBetas ~ normal(0,1);
    ic50 ~ normal(indivMu,indivSigma);
  }
"

assignGroups<-function(x,selector,outId=99999){
  cladeBs<-unique(x[selector])
  notCladeBs<-unique(x[!selector])
  cladeBIds<-structure(c(1:length(cladeBs),rep(outId,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  return(cladeBIds)
}
groupTypes<-sapply(1:max(hiv$group),function(zz)paste(ifelse(hiv[hiv$group==zz,'fluid'][1]=='PL','PL','GE'),ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec')))
groupTypeIds<-structure(1:length(unique(groupTypes)),names=unique(groupTypes))
cladeBIds<-assignGroups(hiv$Pair.ID..,hiv$Subtype=='B')
recipientIds<-assignGroups(hiv$sample[order(hiv$Pair.ID..)],!hiv$donor[order(hiv$Pair.ID..)])
alphaIds<-assignGroups(hiv$sampleSelect[order(hiv$Pair.ID..)],hiv[order(hiv$Pair.ID..),'select']=='A2'&hiv[order(hiv$Pair.ID..),'donor'])
betaIds<-assignGroups(hiv$sampleSelect[order(hiv$Pair.ID..)],hiv[order(hiv$Pair.ID..),'select']=='BE'&hiv[order(hiv$Pair.ID..),'donor'])
recipientAlphaIds<-assignGroups(hiv$sampleSelect[order(hiv$Pair.ID..)],hiv[order(hiv$Pair.ID..),'select']=='A2'&!hiv[order(hiv$Pair.ID..),'donor'])
recipientBetaIds<-assignGroups(hiv$sampleSelect[order(hiv$Pair.ID..)],hiv[order(hiv$Pair.ID..),'select']=='BE'&!hiv[order(hiv$Pair.ID..),'donor'])

fits<-lapply(names(targetCols),function(targetCol){
  #just group by donor or recipient
  #groupTypes<-sapply(1:max(hiv$group),function(zz)ifelse(hiv[hiv$group==zz,'donor'][1],'Don','Rec'))
  #note 99999 is a arbitrarily high number for non clade Bs (should never be called within Stan due to if(cladeB))
  dat<-withAs('xx'=hiv[!is.na(hiv[,targetCol]),],list(
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
    nRecipient=sum(recipientIds<9999),
    recipientIds=recipientIds[xx$sample],
    nGroupTypes=length(unique(groupTypes)),
    groupTypes=groupTypeIds[groupTypes],
    nCladeB=sum(cladeBIds<9999),
    cladeBId=cladeBIds[as.character(xx$Pair.ID..)],
    nAlpha=sum(alphaIds<9999),
    alphaIds=alphaIds[xx$sampleSelect],
    nBeta=sum(betaIds<9999),
    betaIds=betaIds[xx$sampleSelect],
    nRecipientAlpha=sum(recipientAlphaIds<9999),
    recipientAlphaIds=recipientAlphaIds[xx$sampleSelect],
    nRecipientBeta=sum(recipientBetaIds<9999),
    recipientBetaIds=recipientBetaIds[xx$sampleSelect],
    isBeta=as.numeric(xx$select=='BE')
  ))
  fit <- cacheOperation(sprintf('work/stan%s.Rdat',targetCol),stan,model_code = stanCode, data = dat, iter=50000, chains=nThreads,thin=25,control=list(adapt_delta=.999,stepsize=.01))
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
  allPars<-c("metaDonorMu", "metaDonorSd", "metaRecipientMu", "metaRecipientSd", "metaGenitalMu", "metaGenitalSd","metaCladeMu","metaCladeSd","donors", "sigmas", "metaSigmaMu", "metaSigmaSigma", "genitals", "recipients", "clades","metaAlphaMu","metaAlphaSd","metaBetaMu","metaBetaSd","alphas","betas","metaRecipientAlphaMu","metaRecipientAlphaSd","metaRecipientBetaMu","metaRecipientBetaSd","recipientAlphas","recipientBetas")
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
    #pairs(fit,pars=c("metaSigmaMu","metaSigmaSigma", "sigmas"))
    #pairs(fit,pars=c("genitals","metaGenitalMu","metaGenitalSd"))
    #pairs(fit,pars=c("recipients","metaRecipientMu","metaRecipientSd"))
    #pairs(fit,pars=c("clades","metaCladeMu","metaCladeSd"))
    #pairs(fit,pars=c("donors","metaDonorMu","metaDonorSd"))
  dev.off()
  #
  indivRecipBeta<-sims[,grep('recipients\\[[0-9]\\]',colnames(sims))]*sims[,'metaRecipientSd']+sims[,'metaRecipientMu']
  metaBeta<-sims[,'metaRecipientMu']
  indivGenital<-sims[,grep('genitals\\[[0-9]\\]',colnames(sims))]*sims[,'metaGenitalSd']+sims[,'metaGenitalMu']
  metaGenital<-sims[,'metaGenitalMu']
  indivClade<-sims[,grep('clades\\[[0-9]\\]',colnames(sims))]*sims[,'metaCladeSd']+sims[,'metaCladeMu']
  metaClade<-sims[,'metaCladeMu']
  indivAlpha<-sims[,grep('alphas\\[[0-9]\\]',colnames(sims))]*sims[,'metaAlphaSd']+sims[,'metaAlphaMu']
  metaAlpha<-sims[,'metaAlphaMu']
  indivBeta<-sims[,grep('betas\\[[0-9]\\]',colnames(sims))]*sims[,'metaBetaSd']+sims[,'metaBetaMu']
  metaBeta<-sims[,'metaBetaMu']
  indivRecAlpha<-sims[,grep('recipientAlphas\\[[0-9]\\]',colnames(sims))]*sims[,'metaRecipientAlphaSd']+sims[,'metaRecipientAlphaMu']
  metaRecAlpha<-sims[,'metaRecipientAlphaMu']
  indivRecBeta<-sims[,grep('recipientBetas\\[[0-9]\\]',colnames(sims))]*sims[,'metaRecipientBetaSd']+sims[,'metaRecipientBetaMu']
  metaRecBeta<-sims[,'metaRecipientBetaMu']
  metaVarDonor<-log10(sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['PL Don'])])
  metaVarGenital<-log10(sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['GE Don'])])
  metaVarRec<-log10(sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['PL Rec'])])
  indivVarDonor<-log10(sims[,sprintf('sigmas[%d]',which(groupTypes=='PL Don'))]*sims[,sprintf('metaSigmaSigma[%d]',groupTypeIds['PL Don'])]+sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['PL Don'])])
  indivVarGenital<-log10(sims[,sprintf('sigmas[%d]',which(groupTypes=='GE Don'))]*sims[,sprintf('metaSigmaSigma[%d]',groupTypeIds['GE Don'])]+sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['GE Don'])])
  indivVarRec<-log10(sims[,sprintf('sigmas[%d]',which(groupTypes=='PL Rec'))]*sims[,sprintf('metaSigmaSigma[%d]',groupTypeIds['PL Rec'])]+sims[,sprintf('metaSigmaMu[%d]',groupTypeIds['PL Rec'])])
  xlim<-c(-2.7,2.7)
  #
  bins<-seq(xlim[1],xlim[2],.01)
  indivTabs<-apply(indivRecipBeta,2,function(x)table(cut(x,bins))/length(x))
  metaTabs<-table(cut(metaBeta,bins))/length(metaBeta)
  indivGenitalTabs<-apply(indivGenital,2,function(x)table(cut(x,bins))/length(x))
  genitalTabs<-table(cut(metaGenital,bins))/length(metaGenital)
  indivCladeTab<-apply(indivClade,2,function(x)table(cut(x,bins))/length(x))
  cladeTabs<-table(cut(metaClade,bins))/length(metaClade)
  indivAlphaTab<-apply(indivAlpha,2,function(x)table(cut(x,bins))/length(x))
  alphaTabs<-table(cut(metaAlpha,bins))/length(metaClade)
  indivBetaTab<-apply(indivBeta,2,function(x)table(cut(x,bins))/length(x))
  betaTabs<-table(cut(metaBeta,bins))/length(metaClade)
  indivRecAlphaTab<-apply(indivRecAlpha,2,function(x)table(cut(x,bins))/length(x))
  recAlphaTabs<-table(cut(metaRecAlpha,bins))/length(metaClade)
  indivRecBetaTab<-apply(indivRecBeta,2,function(x)table(cut(x,bins))/length(x))
  recBetaTabs<-table(cut(metaRecBeta,bins))/length(metaClade)
  metaVarDonorTab<-table(cut(metaVarDonor,bins))/length(metaVarDonor)
  metaVarGenitalTab<-table(cut(metaVarGenital,bins))/length(metaVarGenital)
  metaVarRecTab<-table(cut(metaVarRec,bins))/length(metaVarRec)
  indivVarDonorTab<-apply(indivVarDonor,2,function(x)table(cut(x,bins))/length(x))
  indivVarGenitalTab<-apply(indivVarGenital,2,function(x)table(cut(x,bins))/length(x))
  indivVarRecTab<-apply(indivVarRec,2,function(x)table(cut(x,bins))/length(x))
  #
  pdf(sprintf('out/bayes/bayes%s.pdf',targetCol),height=9,width=6)
    par(mfrow=c(6,1),las=1,mar=c(3,3.6,1.1,.1))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.7,0),log='x',main='Individual effects',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    meanBin<-(bins[-length(bins)]+bins[-1])/2
    apply(indivTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(indivGenitalTabs,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(indivCladeTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    abline(v=1,lty=2)
    logAxis(1,mgp=c(3,.5,0))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Population effects',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaTabs,0,0),col='#0000FF44',border='#0000FF99')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,genitalTabs,0,0),col='#FF000044',border='#FF000099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,cladeTabs,0,0),col='#00FF0044',border='#00FF0099')
    abline(v=1,lty=2)
    legend('topleft',c('Recipient','Clade B','Genital'),fill=c('#0000FF44','#00FF0044','#FF000044'),border=c('#0000FF99','#00FF0099','#FF000099'),inset=.02)
    logAxis(1,mgp=c(3,.5,0))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Individual selection',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    logAxis(1,mgp=c(3,.5,0))
    abline(v=1,lty=2)
    apply(indivAlphaTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(indivBetaTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(indivRecAlphaTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    apply(indivRecBetaTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FFD70011',border='#FFD70044'))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Population selection',xaxs='i')
    title(xlab=sprintf('Fold increase in %s',gsub('\n',' ',targetCols[targetCol])),mgp=c(1.5,1,0))
    logAxis(1,mgp=c(3,.5,0))
    abline(v=1,lty=2)
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,alphaTabs,0,0),col='#0000FF44',border='#0000FF99')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,betaTabs,0,0),col='#FF000044',border='#FF000099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,recAlphaTabs,0,0),col='#00FF0044',border='#00FF0099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,recBetaTabs,0,0),col='#FFD70099',border='#FFD70099')
    legend('topleft',c('Alpha select','Beta select','Recipient alpha','Recipient beta'),fill=c('#0000FF44','#FF000044','#00FF0044','#FFD70044'),border=c('#0000FF99','#FF000099','#00FF0099','#FFD70099'),inset=.02)
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Individual variance',xaxs='i')
    title(xlab='Standard deviation',mgp=c(1.5,1,0))
    logAxis(1,mgp=c(3,.5,0))
    apply(indivVarDonorTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(indivVarGenitalTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(indivVarRecTab,2,function(xx)polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    plot(1,1,type='n',xlim=10^xlim,ylim=range(indivTabs,metaTabs),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),log="x",main='Population variance',xaxs='i')
    title(xlab='Standard deviation',mgp=c(1.5,1,0))
    logAxis(1,mgp=c(3,.5,0))
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaVarDonorTab,0,0),col='#0000FF44',border='#0000FF99')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaVarGenitalTab,0,0),col='#FF000044',border='#FF000099')
    polygon(10^c(xlim[1],meanBin,xlim[2],xlim[1]),c(0,metaVarRecTab,0,0),col='#00FF0044',border='#00FF0099')
    legend('topleft',c('Donor','Genital','Recipient'),fill=c('#0000FF44','#FF000044','#00FF0044'),border=c('#0000FF99','#FF000099','#00FF0099'),inset=.02)
  dev.off()

}


