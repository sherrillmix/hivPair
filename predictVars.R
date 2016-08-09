mclapply(names(selectVars),function(targetCol){
  message(targetCol)
  modelInput<-as.data.frame(onlyDiff)
  #modelInput$pair<-as.factor(hiv$Pair.ID..)
  modelMatrix<-model.matrix(formula(sprintf('~ %s',paste(colnames(modelInput),collapse='+'))),modelInput)
  selector<-hiv$donor
  target<-log2(hiv[selector,targetCol])-log2(ave(hiv[selector,targetCol],hiv[selector,'Pair.ID..'],FUN=function(x)mean(x,na.rm=TRUE)))
  fitAlpha<-cv.glmnet(modelMatrix[selector,][!is.na(target),],target[!is.na(target)],nfolds=sum(!is.na(target)),grouped=FALSE)
  #fitBeta<-cv.glmnet(modelMatrix[hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.),],log2(hiv$IFNbeta.PD.IC50..ng.ml.[hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.)]),nfolds=sum(hiv$donor&!is.na(hiv$IFNbeta.PD.IC50..ng.ml.)),grouped=FALSE)
  multiFit<-lapply(unique(hiv$Pair.ID..),function(xx){
    pairSelect<-hiv[selector,'Pair.ID..']==xx&!is.na(target)
    cv.glmnet(modelMatrix[selector,][pairSelect,],target[pairSelect],nfolds=sum(pairSelect),grouped=FALSE)
  })
  names(multiFit)<-unique(hiv$Pair.ID..)
  pdf(sprintf('out/lasso/%s.pdf',targetCol))
  plotGlmnet(fitAlpha,markBest1SE=TRUE,main=selectVars[targetCol])
  plotBetas(fitAlpha$glmnet.fit,ylab=selectVars[targetCol],transformFunc=function(x)2^x,labelLambda=fitAlpha$lambda.1se)
  #plotGlmnet(fitBeta,markBest1SE=TRUE,main='IFN beta IC50')
  #plotBetas(fitBeta$glmnet.fit,ylab=expression('IFN beta IC'[50]%*%' '),transformFunc=function(x)2^x,labelLambda=10^-.3)
  coefs<-coef(fitAlpha)
  coefs<-coefs[coefs[,1]!=0,]
  coefs<-coefs[names(coefs)!="(Intercept)"]
  cols<-rainbow.lab(length(unique(hiv$Pair.ID..)),alpha=.7)
  names(cols)<-unique(hiv$Pair.ID..)
  if(length(coefs)>0){
    for(ii in names(coefs)){
      vpPlot(modelInput[selector,sub('[A-Z]$','',ii)][!is.na(target)],target[!is.na(target)],col=NA,bg=cols[as.character(hiv$Pair.ID..[selector])],pch=21,main=sprintf('%s: %0.3f',ii,coefs[ii]),ylab=sprintf('Mean adjusted %s (log)',targetCol))
      legend('topleft',names(cols),pt.bg=cols,pch=21,col=NA)
    }
  }
  for(ii in names(multiFit)){
    message(ii)
    plotGlmnet(multiFit[[ii]],markBest1SE=TRUE,main=sprintf('Pair %s %s',ii,sub('\n',' ',selectVars[targetCol])))
    plotBetas(multiFit[[ii]]$glmnet.fit,ylab=selectVars[targetCol],transformFunc=function(x)2^x,labelLambda=multiFit[[ii]]$lambda.1se)
  }
  dev.off()
},mc.cores=3)

