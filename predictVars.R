library(glmnet)
library(glmnetPlotR)
if(!exists('hiv'))source('readData.R')

out<-mclapply(names(targetCols),function(targetCol){
  message(targetCol)
  modelInput<-as.data.frame(onlyDiffAA)
  #modelInput$pair<-as.factor(hiv$Pair.ID..)
  modelMatrix<-model.matrix(formula(sprintf('~ %s',paste(colnames(modelInput),collapse='+'))),modelInput)
  selector<-hiv$donor&!is.na(hiv[,targetCol])
  unadjustTarget<-log2(hiv[selector,targetCol])
  target<-unadjustTarget-log2(ave(hiv[selector,targetCol],hiv[selector,'Pair.ID..'],FUN=function(x)mean(x,na.rm=TRUE)))
  fitAlpha<-cv.glmnet(modelMatrix[selector,],target,nfolds=length(target),grouped=FALSE)
  fitAlpha2<-cv.glmnet(modelMatrix[selector,],unadjustTarget,nfolds=length(target),grouped=FALSE)
  #fitBeta<-cv.glmnet(modelMatrix[hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.),],log2(hiv$IFNbeta.PD.IC50..ng.ml.[hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.)]),nfolds=sum(hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.)),grouped=FALSE)
  multiFit<-lapply(unique(hiv$Pair.ID..),function(xx){
    pairSelect<-hiv[selector,'Pair.ID..']==xx
    if(sum(pairSelect)<5||length(unique(target[pairSelect]))==1)return(NULL)
    cv.glmnet(modelMatrix[selector,][pairSelect,],target[pairSelect],nfolds=sum(pairSelect),grouped=FALSE)
  })
  names(multiFit)<-unique(hiv$Pair.ID..)
  pdf(sprintf('out/lasso/%s.pdf',targetCol))
  plotGlmnet(fitAlpha2,markBest1SE=TRUE,main=targetCols[targetCol])
  plotBetas(fitAlpha2$glmnet.fit,ylab=targetCols[targetCol],transformFunc=function(x)2^x,labelLambda=fitAlpha$lambda.1se)
  coefs<-coef(fitAlpha2)
  coefs<-coefs[coefs[,1]!=0,]
  coefs<-coefs[names(coefs)!="(Intercept)"]
  cols<-rainbow.lab(length(unique(hiv$Pair.ID..)),alpha=.7)
  names(cols)<-unique(hiv$Pair.ID..)
  if(length(coefs)>0){
    for(ii in names(coefs)){
      vpPlot(modelInput[selector,sub('[A-Z]$','',ii)],target,col=NA,bg=cols[as.character(hiv$Pair.ID..[selector])],pch=21,main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('Mean adjusted %s (log)',targetCol))
      legend('topleft',names(cols),pt.bg=cols,pch=21,col=NA)
    }
  }
  plotGlmnet(fitAlpha,markBest1SE=TRUE,main=targetCols[targetCol])
  plotBetas(fitAlpha$glmnet.fit,ylab=targetCols[targetCol],transformFunc=function(x)2^x,labelLambda=fitAlpha$lambda.1se)
  #plotGlmnet(fitBeta,markBest1SE=TRUE,main='IFN beta IC50')
  #plotBetas(fitBeta$glmnet.fit,ylab=expression('IFN beta IC'[50]%*%' '),transformFunc=function(x)2^x,labelLambda=10^-.3)
  coefs<-coef(fitAlpha)
  coefs<-coefs[coefs[,1]!=0,]
  coefs<-coefs[names(coefs)!="(Intercept)"]
  cols<-rainbow.lab(length(unique(hiv$Pair.ID..)),alpha=.7)
  names(cols)<-unique(hiv$Pair.ID..)
  if(length(coefs)>0){
    for(ii in names(coefs)){
      vpPlot(modelInput[selector,sub('[A-Z]$','',ii)],target,col=NA,bg=cols[as.character(hiv$Pair.ID..[selector])],pch=21,main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('Mean adjusted %s (log)',targetCol))
      legend('topleft',names(cols),pt.bg=cols,pch=21,col=NA)
    }
  }
  for(ii in names(multiFit)){
    if(is.null(multiFit[[ii]]))next()
    message(ii)
    plotGlmnet(multiFit[[ii]],markBest1SE=TRUE,main=sprintf('Pair %s %s',ii,sub('\n',' ',targetCols[targetCol])))
    plotBetas(multiFit[[ii]]$glmnet.fit,ylab=targetCols[targetCol],transformFunc=function(x)2^x,labelLambda=multiFit[[ii]]$lambda.1se)
    coefs<-coef(multiFit[[ii]])
    coefs<-coefs[coefs[,1]!=0,]
    coefs<-coefs[names(coefs)!="(Intercept)"]
    if(length(coefs)>0){
      for(ii in names(coefs)){
        vpPlot(modelInput[selector,sub('[A-Z]$','',ii)][!is.na(unadjustTarget)],unadjustTarget[!is.na(unadjustTarget)],col=c(NA,'black')[(hiv$Pair.ID..[selector]==ii)+1],bg=cols[as.character(hiv$Pair.ID..[selector])],pch=21,main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('%s (log)',targetCol))
        legend('topleft',names(cols),pt.bg=cols,pch=21,col=NA)
      }
    }
  }
  dev.off()
},mc.cores=4)

