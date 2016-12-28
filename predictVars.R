library(glmnet)
library(glmnetPlotR)
if(!exists('onlyDiffAA'))source('parseSeqs.R')

modelInput<-as.data.frame(cbind(onlyDiffAA,onlyDiff))
#modelInput<-as.data.frame(cbind(onlyDiffAA))
#modelInput$pair<-as.factor(hiv$Pair.ID)
modelMatrix<-model.matrix(formula(sprintf('~ %s',paste(colnames(modelInput),collapse='+'))),modelInput)

lassoCols<-c(
  targetCols,
  'Autologous.IC50'='Autologous IC50',
  'Bnaber.IC50'='Bnaber IC50'
)
targetColTransform['Autologous.IC50']<-'log'
targetColTransform['Bnaber.IC50']<-'log'

#arbitrarily take first recipient virus
hiv$firstRecipient<-FALSE
hiv[!hiv$donor&hiv$select=='UT','firstRecipient']<-ave(hiv[!hiv$donor&hiv$select=='UT','firstRecipient'],hiv[!hiv$donor&hiv$select=='UT','baseName'],FUN=function(x)c(TRUE,rep(FALSE,length(x)-1)))

aaFits<-mclapply(names(lassoCols),function(targetCol){
  message(targetCol)
  thisTransform<-targetColTransform[targetCol]
  if(thisTransform=='log'){
    transformFunc<-log10
  }else if(thisTransform=='logit'){
    transformFunc<-function(x)logit(x/100)
  } else if(thisTransform=='identity'){
    transformFunc<-function(x)x
  }else{
    stop(simpleError('Unknown tranform'))
  }
  #selector<-hiv$donor&!is.na(hiv[,targetCol])
  selector<-!is.na(hiv[,targetCol])&(hiv$donor|hiv$firstRecipient)
  #condense columns with no differences
  collapsedCols<-ave(colnames(modelMatrix[selector,]),apply(modelMatrix,2,paste,collapse=''),FUN=function(xx)paste(xx,collapse=','))
  collapsedModel<-modelMatrix
  colnames(collapsedModel)<-collapsedCols
  collapsedModel<-collapsedModel[,!duplicated(collapsedCols)]
  #raw values 
  unadjustTarget<-transformFunc(hiv[selector,targetCol])
  target<-unadjustTarget-transformFunc(ave(hiv[selector,targetCol],hiv[selector,'Pair.ID'],FUN=function(x)mean(x,na.rm=TRUE)))
  fitAlpha<-cv.glmnet(collapsedModel[selector,],target,nfolds=length(target),grouped=FALSE)
  fitAlpha2<-cv.glmnet(collapsedModel[selector,],unadjustTarget,nfolds=length(target),grouped=FALSE)
  multiFit<-lapply(unique(hiv$Pair.ID),function(xx){
    pairSelect<-hiv[selector,'Pair.ID']==xx
    if(sum(pairSelect)<5||length(unique(target[pairSelect]))==1)return(NULL)
    cv.glmnet(collapsedModel[selector,][pairSelect,],target[pairSelect],nfolds=sum(pairSelect),grouped=FALSE)
  })
  names(multiFit)<-unique(hiv$Pair.ID)
  return(list('fitAlpha'=fitAlpha,'fitAlpha2'=fitAlpha2,'multiFit'=multiFit,'target'=target,'unadjustTarget'=unadjustTarget,'selector'=selector))
},mc.cores=5)
names(aaFits)<-names(lassoCols)

minBeta<-.01
for(targetCol in names(lassoCols)){
  thisTransform<-targetColTransform[targetCol]
  thisFit<-aaFits[[targetCol]]
  fitAlpha<-thisFit[['fitAlpha']];fitAlpha2<-thisFit[['fitAlpha2']];multiFit<-thisFit[['multiFit']];target<-thisFit[['target']];unadjustTarget<-thisFit[['unadjustTarget']];selector<-thisFit[['selector']]
  pdf(sprintf('out/lasso/%s.pdf',targetCol))
  plotGlmnet(fitAlpha2,markBest1SE=TRUE,main=lassoCols[targetCol])
  plotBetas(fitAlpha2$glmnet.fit,ylab=lassoCols[targetCol],transformFunc=function(x)2^x,labelLambda=fitAlpha$lambda.1se,minBeta=minBeta)
  coefs<-coef(fitAlpha2)
  coefs<-coefs[abs(coefs[,1])>minBeta,]
  coefs<-coefs[names(coefs)!="(Intercept)"]
  coefs<-coefs[order(abs(coefs),decreasing=TRUE)]
  cols<-rainbow.lab(length(unique(hiv$Pair.ID)),alpha=.7)
  names(cols)<-unique(hiv$Pair.ID)
  if(length(coefs)>0){
    for(ii in names(coefs)){
      par(mar=c(4,4,1.5,.5))
      vpPlot(apply(modelInput[selector,sub('[A-Z]$','',strsplit(ii,',')[[1]]),drop=FALSE],1,paste,collapse=''),unadjustTarget,col=NA,bg=cols[as.character(hiv$Pair.ID[selector])],pch=(21:22)[(!hiv$donor[selector])+1],main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('%s (%s)',targetCol,thisTransform),las=2)
      legend('topleft',c(names(cols),'D','R'),pt.bg=c(cols,'#00000066','#00000066'),pch=rep(21:22,c(length(cols)+1,1)),col=NA)
    }
  }
  plotGlmnet(fitAlpha,markBest1SE=TRUE,main=lassoCols[targetCol])
  plotBetas(fitAlpha$glmnet.fit,ylab=sprintf('Mean adjusted %s',lassoCols[targetCol]),transformFunc=function(x)2^x,labelLambda=fitAlpha$lambda.1se,minBeta=minBeta)
  #plotGlmnet(fitBeta,markBest1SE=TRUE,main='IFN beta IC50')
  #plotBetas(fitBeta$glmnet.fit,ylab=expression('IFN beta IC'[50]%*%' '),transformFunc=function(x)2^x,labelLambda=10^-.3)
  coefs<-coef(fitAlpha)
  coefs<-coefs[abs(coefs[,1])>minBeta,]
  coefs<-coefs[names(coefs)!="(Intercept)"]
  coefs<-coefs[order(abs(coefs),decreasing=TRUE)]
  cols<-rainbow.lab(length(unique(hiv$Pair.ID)),alpha=.7)
  names(cols)<-unique(hiv$Pair.ID)
  if(length(coefs)>0){
    for(ii in names(coefs)){
      par(mar=c(4,4,1.5,.5))
      vpPlot(apply(modelInput[selector,sub('[A-Z]$','',strsplit(ii,',')[[1]]),drop=FALSE],1,paste,collapse=''),target,col=NA,bg=cols[as.character(hiv$Pair.ID[selector])],pch=(21:22)[(!hiv$donor[selector])+1],main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('Mean adjusted %s (%s)',targetCol,thisTransform),las=2)
      legend('topleft',c(names(cols),'D','R'),pt.bg=c(cols,'#00000066','#00000066'),pch=rep(21:22,c(length(cols)+1,1)),col=NA)
    }
  }
  for(pair in sort(names(multiFit))){
    if(is.null(multiFit[[pair]]))next()
    message(pair)
    plotGlmnet(multiFit[[pair]],markBest1SE=TRUE,main=sprintf('Pair %s %s',pair,sub('\n',' ',lassoCols[targetCol])))
    plotBetas(multiFit[[pair]]$glmnet.fit,ylab=lassoCols[targetCol],transformFunc=function(x)10^x,labelLambda=multiFit[[pair]]$lambda.1se,minBeta=minBeta)
    coefs<-coef(multiFit[[pair]])
    coefs<-coefs[abs(coefs[,1])>minBeta,]
    coefs<-coefs[names(coefs)!="(Intercept)"]
    coefs<-coefs[order(abs(coefs),decreasing=TRUE)]
    if(length(coefs)>0){
      for(ii in names(coefs)){
        par(mar=c(4,4,1.5,.5))
        vpPlot(apply(modelInput[selector,sub('[A-Z]$','',strsplit(ii,',')[[1]]),drop=FALSE],1,paste,collapse=''),unadjustTarget,col=c(NA,'black')[(hiv$Pair.ID[selector]==pair)+1],bg=cols[as.character(hiv$Pair.ID[selector])],pch=(21:22)[(!hiv$donor[selector])+1],main=sprintf('Pair %s %s: %0.3f',pair,ii,coefs[ii]),ylab=sprintf('%s (%s)',targetCol,thisTransform),las=2)
        legend('topleft',c(names(cols),'D','R'),pt.bg=c(cols,'#00000066','#00000066'),pch=rep(21:22,c(length(cols)+1,1)),col=c(c(NA,'black')[(names(cols)==pair)+1],NA,NA))
      }
    }
  }
  dev.off()
}

