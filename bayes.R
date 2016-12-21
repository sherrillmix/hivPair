library('rstan')
library('vioplot') #for violin plot in check figure
library('png') #for raster inside pdf
library('parallel')

#set up parallel options
nThreads<-5 #used 50 for final figures
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if(!dir.exists(file.path('out','bayes')))dir.create(file.path('out','bayes'))

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
    int<lower=1> genitalIds[N];
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
      if(isGenital[ii])indivMu[ii]=indivMu[ii]+metaGenitalMu+genitals[genitalIds[ii]]*metaGenitalSd;
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

assignGroups<-function(x,selector=rep(TRUE,length(x)),outId=99999){
  cladeBs<-unique(x[selector])
  notCladeBs<-unique(x[!selector])
  cladeBIds<-structure(c(1:length(cladeBs),rep(outId,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  return(cladeBIds)
}
cladeBIds<-assignGroups(hiv[,'Pair ID'],hiv$Subtype=='B')
recipientIds<-assignGroups(hiv$sample[order(hiv[,'Pair ID'])],!hiv$isDonor[order(hiv[,'Pair ID'])])
alphaIds<-assignGroups(hiv$sampleSelect[order(hiv[,'Pair ID'])],hiv[order(hiv[,'Pair ID']),'Selection']=='A2'&hiv[order(hiv[,'Pair ID']),'isDonor'])
betaIds<-assignGroups(hiv$sampleSelect[order(hiv[,'Pair ID'])],hiv[order(hiv[,'Pair ID']),'Selection']=='BE'&hiv[order(hiv[,'Pair ID']),'isDonor'])
recipientAlphaIds<-assignGroups(hiv$sampleSelect[order(hiv[,'Pair ID'])],hiv[order(hiv[,'Pair ID']),'Selection']=='A2'&!hiv[order(hiv[,'Pair ID']),'isDonor'])
recipientBetaIds<-assignGroups(hiv$sampleSelect[order(hiv[,'Pair ID'])],hiv[order(hiv[,'Pair ID']),'Selection']=='BE'&!hiv[order(hiv[,'Pair ID']),'isDonor'])
genitalIds<-assignGroups(hiv$sampleFluid[order(hiv[,'Pair ID'])],hiv$isGenital[order(hiv[,'Pair ID'])])

#stick recipient beta/alpha treated with untreated recipient variance
hiv$group<-paste(hiv$sample,ifelse(hiv$isGenital,'GE','PL'),ifelse(hiv$isDonor,hiv$Selection,'UT'))
hiv$groupType<-paste(ifelse(hiv$isDonor,'Donor','Recipient'),ifelse(hiv$isGenital,'GE','PL'),ifelse(hiv$isDonor,hiv$Selection,'UT'))
groupIds<-assignGroups(hiv$group[order(hiv[,'Pair ID'],hiv$sampleFluid)],rep(TRUE,nrow(hiv)))
groupTypeIds<-assignGroups(hiv$groupType[order(hiv[,'Pair ID'],hiv$sampleFluid)],rep(TRUE,nrow(hiv)))
groupTypes<-sapply(names(groupIds),function(x)hiv[hiv$group==x,'groupType'][1])

if(!exists('fits')){
  fits<-lapply(targetCols,function(targetCol){
    thisTransform<-targetColTransform[targetCol]
    thisCensorDown<-targetColCensorDown[[targetCol]]
    if(thisTransform=='log'){
      transformFunc<-log10 
    }else if(thisTransform=='logit'){
      transformFunc<-function(x)logit(x/100)
    } else if(thisTransform=='identity'){
      transformFunc<-function(x)x
    }else{
      stop(simpleError('Unknown tranform'))
    }
    xx<-hiv[!is.na(hiv[,targetCol]),]
    dat<-list(
      ic50=transformFunc(xx[,targetCol]),
      N=nrow(xx),
      isCladeB=as.integer(xx$Subtype=='B'),
      isGenital=as.integer(xx$Fluid!='PL'),
      nGenital=sum(genitalIds<9999),
      isRecipient=as.integer(!xx$isDonor),
      nGroup=max(groupIds),
      groupIds=groupIds[xx$group],
      nPair=max(xx[,'Pair ID']),
      pairIds=xx[,'Pair ID'],
      nRecipient=sum(recipientIds<9999),
      genitalIds=genitalIds[xx$sampleFluid],
      recipientIds=recipientIds[xx$sample],
      nGroupTypes=length(unique(xx$groupType)),
      groupTypes=groupTypeIds[groupTypes],
      nCladeB=sum(cladeBIds<9999),
      cladeBId=cladeBIds[as.character(xx[,'Pair ID'])],
      nAlpha=sum(alphaIds<9999),
      alphaIds=alphaIds[xx$sampleSelect],
      nBeta=sum(betaIds<9999),
      betaIds=betaIds[xx$sampleSelect],
      nRecipientAlpha=sum(recipientAlphaIds<9999),
      recipientAlphaIds=recipientAlphaIds[xx$sampleSelect],
      nRecipientBeta=sum(recipientBetaIds<9999),
      recipientBetaIds=recipientBetaIds[xx$sampleSelect],
      isBeta=as.numeric(xx$Selection=='BE') #don't actually use this but keeping in for now for caching
    )
    #R makes local copy so we don't have to worry about overwrite global stanCode
    if(thisTransform=='identity'){
      stanCode<-gsub('~ *gamma\\([0-9,]+\\)','~ gamma(2,.01)',stanCode)
    }
    if(!is.null(thisCensorDown)){
      cutVal<-transformFunc(thisCensorDown[!is.na(hiv[,targetCol])])
      #multiply by a small amount to avoid any floating point weirdness with comparison 
      stanCode<-sub(
        'ic50 ~ normal\\(indivMu,indivSigma\\);',
        sprintf('for (ii in 1:N){\nif(ic50[ii]<= censorDown[ii]*(1+1e-6))target+=normal_lcdf(censorDown[ii]|indivMu[ii],indivSigma[ii]);\nelse ic50[ii] ~ normal(indivMu[ii],indivSigma[ii]); \n}'),
        stanCode
      )
      stanCode<-sub('real ic50\\[N\\];','real ic50[N];\nreal censorDown[N];',stanCode)
      dat<-c(dat,list('censorDown'=cutVal))
    }
    #used iter=100000 and thin=25 for final figures
    fit <- stan(model_code = stanCode, data = dat, iter=4000, chains=nThreads,thin=10,control=list(adapt_delta=.999,stepsize=.01))
    return(list('fit'=fit,'dat'=dat,stan=stanCode))
  })
  names(fits)<-targetCols
}

#convert the N(0,1) dists to real distributions based on meta mu and sd
convertCols<-function(cols,means,sds,sims){
  out<-mapply(function(col,mean,sd){
    if(length(col)==1)col<-grep(col,colnames(sims))
    sims[,col]*sims[,sd]+sims[,mean]
  },cols,means,sds,SIMPLIFY=FALSE) 
  return(do.call(cbind,out))
}

cachedTabs<-lapply(targetCols,function(targetCol){
  if(targetColTransform[targetCol]=='identity') xlim<-c(-420,530)
  else if(targetCol=='vres')xlim<-c(-1.2,3.5)
  else xlim<-c(-1.2,2.5)
  bins<-seq(xlim[1],xlim[2],length.out=200)
  message(targetCol)
  dat<-fits[[targetCol]][['dat']]
  fit<-fits[[targetCol]][['fit']]
  sims<-as.array(fit)
  dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
  colnames(sims)<-dimnames(as.array(fit))[[3]]
  #
  converted<-convertCols(
    list('recipients\\[[0-9]+\\]','genitals\\[[0-9]+\\]','clades\\[[0-9]+\\]','alphas\\[[0-9]+\\]','betas\\[[0-9]+\\]','recipientAlphas\\[[0-9]\\]','recipientBetas\\[[0-9]\\]',sprintf('sigmas[%d]',1:max(dat$groupIds))),
    list('metaRecipientMu','metaGenitalMu','metaCladeMu','metaAlphaMu','metaBetaMu','metaRecipientAlphaMu','metaRecipientBetaMu',sprintf('metaSigmaMu[%d]',dat$groupTypes)),
    list('metaRecipientSd','metaGenitalSd','metaCladeSd','metaAlphaSd','metaBetaSd','metaRecipientAlphaSd','metaRecipientBetaSd',sprintf('metaSigmaSigma[%d]',dat$groupTypes)),
    sims
  )
  converted<-cbind(converted,sims[,c('metaDonorMu','metaRecipientMu','metaGenitalMu','metaCladeMu','metaAlphaMu','metaBetaMu','metaRecipientAlphaMu','metaRecipientBetaMu',sprintf('metaSigmaMu[%d]',1:max(dat$groupTypes)))])
  converted[,grep('^(sigmas|metaSigmaMu)\\[[0-9]+\\]',colnames(converted))]<-log10(converted[,grep('^(sigmas|metaSigmaMu)\\[[0-9]+\\]',colnames(converted))])
  donorFolds<-converted[,colnames(converted)[-grep('\\[',colnames(converted))]]/converted[,'metaDonorMu']
  colnames(donorFolds)<-sprintf('fold_%s',colnames(donorFolds))
  converted<-cbind(converted,donorFolds)
  stats<-apply(converted,2,function(x)c('mean'=mean(x),quantile(x,c(.025,.05,.95,.975)),'gt0'=mean(x>0),'lt0'=mean(x<0),n=length(x)))
  stats2<-c(
    'p(beta>alpha)'=mean(converted[,'metaBetaMu']>converted[,'metaAlphaMu']),
    'p(|recipient-beta|<|recipient-alpha|)'=mean(abs(converted[,'metaRecipientMu']-converted[,'metaBetaMu'])<abs(converted[,'metaRecipientMu']-converted[,'metaAlphaMu'])),
    'p(recipient>alpha)'=mean(converted[,'metaRecipientMu']>converted[,'metaAlphaMu']),
    'p(recipient>beta)'=mean(converted[,'metaRecipientMu']>converted[,'metaBetaMu']),
    'p(beta<alpha)'=mean(converted[,'metaBetaMu']<converted[,'metaAlphaMu']),
    'p(|recipient-beta|>|recipient-alpha|)'=mean(abs(converted[,'metaRecipientMu']-converted[,'metaBetaMu'])>abs(converted[,'metaRecipientMu']-converted[,'metaAlphaMu'])),
    'p(recipient<alpha)'=mean(converted[,'metaRecipientMu']<converted[,'metaAlphaMu']),
    'p(recipient<beta)'=mean(converted[,'metaRecipientMu']<converted[,'metaBetaMu'])
  )
  tabbed<-apply(converted,2,function(x)table(cut(x,bins))/length(x))
  return(list('tabs'=tabbed,'stats'=stats,'stats2'=stats2,'bins'=bins))
})
names(cachedTabs)<-targetCols

allPars<-c("metaDonorMu", "metaDonorSd", "metaRecipientMu", "metaRecipientSd", "metaGenitalMu", "metaGenitalSd","metaCladeMu","metaCladeSd","donors", "sigmas", "metaSigmaMu", "metaSigmaSigma", "genitals", "recipients", "clades","metaAlphaMu","metaAlphaSd","metaBetaMu","metaBetaSd","alphas","betas","metaRecipientAlphaMu","metaRecipientAlphaSd","metaRecipientBetaMu","metaRecipientBetaSd","recipientAlphas","recipientBetas")

#generate figures and stats
for(targetCol in targetCols){
  message(targetCol)
  #
  if(targetColTransform[targetCol]=='identity'){
    transform<-function(x)x
    logX<-''
    xlab<-sprintf('Increase in %s',gsub('\n',' ',targetCol))
  }
  else{
    transform<-function(x)10^x
    logX<-'x'
    xlab<-sprintf('Fold increase in %s',gsub('\n',' ',targetCol))
    if(targetColTransform[targetCol]=='logit')xlab<-sprintf('%s odds',xlab)
  }
  #
  fit<-fits[[targetCol]][['fit']]
  dat<-fits[[targetCol]][['dat']]
  tabs<-cachedTabs[[targetCol]][['tabs']]
  stats<-cachedTabs[[targetCol]][['stats']]
  stats2<-cachedTabs[[targetCol]][['stats2']]
  bins<-cachedTabs[[targetCol]][['bins']]
  xlim<-range(bins)
  #
  recipientCols<-grep('recipients\\[[0-9]+\\]',colnames(tabs))
  genitalCols<-grep('genitals\\[[0-9]+\\]',colnames(tabs))
  cladeCols<-grep('clades\\[[0-9]+\\]',colnames(tabs))
  alphaCols<-grep('alphas\\[[0-9]+\\]',colnames(tabs))
  betaCols<-grep('betas\\[[0-9]+\\]',colnames(tabs))
  recAlphaCols<-grep('recipientAlphas\\[[0-9]\\]',colnames(tabs))
  recBetaCols<-grep('recipientBetas\\[[0-9]\\]',colnames(tabs))
  metaVarCols<-structure(sprintf('metaSigmaMu[%d]',groupTypeIds),names=names(groupTypeIds))
  varCols<-sprintf('sigmas[%d]',1:max(dat$groupIds))
  varColsByGroup<-tapply(varCols,groupTypes,c)
  #
  outStatCols<-c(
    'Recipient fold change'='metaRecipientMu',
    'Genital fold change'='metaGenitalMu',
    'Clade B fold change'='metaCladeMu',
    'Alpha selection fold change'='metaAlphaMu',
    'Beta selection fold change'='metaBetaMu',
    'Recipient alpha selection fold change'='metaRecipientAlphaMu',
    'Recipient beta selection fold change'='metaRecipientBetaMu'
  )
  if(targetColTransform[targetCol]=='identity'){
    folds<-sprintf('fold_%s',outStatCols)
    names(outStatCols)<-sub('fold ','',names(outStatCols))
    outStatCols<-c(outStatCols,folds)
  }
  outStats<-as.data.frame(t(stats[,outStatCols]))
  outStats$mean<-transform(outStats$mean)
  outStats$'95% CrI'<-sprintf('%s-%s',sapply(signif(transform(outStats[,'2.5%']),digits=3),formatC,digits=3,format='fg',flag='#'),sapply(signif(transform(outStats[,'97.5%']),digits=3),formatC,digits=3,format='fg',flag='#'))
  outStats$'90% CrI'<-sprintf('%s-%s',sapply(signif(transform(outStats[,'5%']),digits=3),formatC,digits=3,format='fg',flag='#'),sapply(signif(transform(outStats[,'95%']),digits=3),formatC,digits=3,format='fg',flag='#'))
  rownames(outStats)<-names(outStatCols)
  outStats$'p(effect<=1)'<-format(1-outStats$gt0,digits=3)
  outStats[outStats$gt0==1,'p(effect<=1)']<-sprintf("<%s",format(1/outStats[outStats$gt0==1,'n'],digits=1,scientific=FALSE))
  output<-outStats[,c('mean','95% CrI','90% CrI','p(effect<=1)')]
  if(targetColTransform[targetCol]=='identity')colnames(output)[colnames(output)=='p(effect<=1)']<-'p(effect<=0)'
  write.csv(output,file.path('out','bayes',sprintf('stats_%s.csv',sub('/','_',targetCol))))
  outStats2<-data.frame('probability'=stats2)
  outStats2[outStats2$probability==0,'probability']<-sprintf("<%s",format(1/outStats[1,'n'],digits=1,scientific=FALSE))
  write.csv(outStats2,file.path('out','bayes',sprintf('stats2_%s.csv',sub('/','_',targetCol))))
  #
  tmp<-cladeBIds
  names(tmp)<-sprintf('Clade %s',names(cladeBIds))
  #just using genitals directly since first 3 were genitals but really should use an id for genital instead
  tmp2<-structure(1:100,names=sprintf('Genital %d',1:100))
  converts<-list('recipients'=recipientIds,'clades'=tmp,'alphas'=alphaIds,'betas'=betaIds,'recipientAlphas'=recipientAlphaIds,'recipientBetas'=recipientBetaIds,'genitals'=tmp2)
  #
  pdf(file.path('out','bayes',sprintf('bayes%s.pdf',sub('/','_',targetCol))),height=9,width=6)
    par(mfrow=c(6,1),las=1,mar=c(3,3.6,1.1,.1))
    ylims<-c(0,max(tabs[,c(colnames(tabs)[c(recipientCols,cladeCols,genitalCols,alphaCols,betaCols)],'metaRecipientMu','metaGenitalMu','metaCladeMu','metaAlphaMu','metaBetaMu')]))
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.7,0),log=logX,main='Individual effects',xaxs='i')
    title(xlab=xlab,mgp=c(1.5,1,0))
    meanBin<-(bins[-length(bins)]+bins[-1])/2
    apply(tabs[,recipientCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(tabs[,genitalCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(tabs[,cladeCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    abline(v=1,lty=2)
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.8,0),log=logX,main='Population effects',xaxs='i')
    title(xlab=xlab,mgp=c(1.5,1,0))
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaRecipientMu'],0,0),col='#0000FF44',border='#0000FF99')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaGenitalMu'],0,0),col='#FF000044',border='#FF000099')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaCladeMu'],0,0),col='#00FF0044',border='#00FF0099')
    abline(v=1,lty=2)
    legend('topleft',c('Recipient','Clade B','Genital'),fill=c('#0000FF44','#00FF0044','#FF000044'),border=c('#0000FF99','#00FF0099','#FF000099'),inset=.02)
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.8,0),log=logX,main='Individual selection',xaxs='i')
    title(xlab=xlab,mgp=c(1.5,1,0))
    abline(v=1,lty=2)
    apply(tabs[,alphaCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(tabs[,betaCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(tabs[,recAlphaCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    apply(tabs[,recBetaCols],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#FFD70011',border='#FFD70044'))
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.8,0),log=logX,main='Population selection',xaxs='i')
    title(xlab=xlab,mgp=c(1.5,1,0))
    abline(v=1,lty=2)
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaRecipientMu'],0,0),col='#00000022',border='#00000099',lty=2)
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaAlphaMu'],0,0),col='#0000FF44',border='#0000FF99')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaBetaMu'],0,0),col='#FF000044',border='#FF000099')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaRecipientAlphaMu'],0,0),col='#00FF0044',border='#00FF0099')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaRecipientBetaMu'],0,0),col='#FFD70044',border='#FFD70099')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,'metaGenitalMu'],0,0),col='#CC33FF44',border='#CC33FF99')
    legend('topleft',c('Alpha select','Beta select','Recipient alpha','Recipient beta','Untreated recipient','Genital'),fill=c('#0000FF44','#FF000044','#00FF0044','#FFD70044','#00000022','#CC33FF44'),border=c('#0000FF99','#FF000099','#00FF0099','#FFD70099','#00000099','#CC33FF99'),inset=.02)
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.8,0),log=logX,main='Individual variance',xaxs='i')
    title(xlab='Standard deviation',mgp=c(1.5,1,0))
    apply(tabs[,varColsByGroup[['Donor PL UT']]],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#0000FF11',border='#0000FF44'))
    apply(tabs[,varColsByGroup[['Donor GE UT']]],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#FF000011',border='#FF000044'))
    apply(tabs[,varColsByGroup[['Recipient PL UT']]],2,function(xx)polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,xx,0,0),col='#00FF0011',border='#00FF0044'))
    plot(1,1,type='n',xlim=transform(xlim),ylim=ylims,xlab='',ylab='Posterior probability',mgp=c(2.7,.8,0),log=logX,main='Population variance',xaxs='i')
    title(xlab='Standard deviation',mgp=c(1.5,1,0))
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,metaVarCols['Donor PL UT']],0,0),col='#0000FF44',border='#0000FF99')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,metaVarCols['Donor GE UT']],0,0),col='#FF000044',border='#FF000099')
    polygon(transform(c(xlim[1],meanBin,xlim[2],xlim[1])),c(0,tabs[,metaVarCols['Recipient PL UT']],0,0),col='#00FF0044',border='#00FF0099')
    legend('topleft',c('Donor','Genital','Recipient'),fill=c('#0000FF44','#FF000044','#00FF0044'),border=c('#0000FF99','#FF000099','#00FF0099'),inset=.02)
  dev.off()
}

#generate figures to check convergence and posterior predictives
check<-mclapply(targetCols,function(targetCol){
  message(targetCol)
  fit<-fits[[targetCol]][['fit']]
  dat<-fits[[targetCol]][['dat']]
  #
  sims<-as.array(fit)
  dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
  colnames(sims)<-dimnames(as.array(fit))[[3]]
  #
  indivMuCols<-sprintf('indivMu[%d]',1:dat[['N']])
  indivSdCols<-sprintf('indivSigma[%d]',1:dat[['N']])
  simFits<-mapply(function(mu,sigma)rnorm(nrow(sims),sims[,mu],sims[,sigma]),indivMuCols,indivSdCols,SIMPLIFY=FALSE)
  #
  tmp<-tempfile()
  png(file=tmp,width=2000,height=2000,res=100)
    print(traceplot(fit,pars=allPars))
  dev.off()
  tmp2<-tempfile()
  png(file=tmp2,width=2000,height=2000,res=100)
    par(mar=c(4,5,2,0))
    plot(1,1,type='n',xlim=c(1,dat[['N']])+c(-1,1),ylim=range(unlist(simFits)),xlab='Virus ID',ylab=sprintf('Log10 %s',targetCol),xaxs='i',cex.axis=1.5,cex.lab=2,main='Posterior predictive distributions',cex.main=2)
    do.call(vioplot,c('x'=list(simFits[[1]]),simFits[-1],'add'=list(TRUE),'colMed'=list(NA),'col'=list('#00000033')))
    points(1:length(dat[['ic50']]),dat[['ic50']],col='red',pch='+',lwd=2)
  dev.off()
  pdf(file.path('out','bayes',sprintf('fit%s.pdf',sub('/','_',targetCol))),width=20,height=20)
    print(plot(fit,pars=allPars)+theme_minimal(base_family="Helvetica"))
    #use rasters for file size
    par(mai=c(0,0,0,0))
    plot(c(0,1),c(0,1),type="n")
    rasterImage(readPNG(tmp),0,0,1,1)
    plot(c(0,1),c(0,1),type="n")
    rasterImage(readPNG(tmp2),0,0,1,1)
  dev.off()
  file.remove(tmp)
  file.remove(tmp2)
},mc.cores=3)
