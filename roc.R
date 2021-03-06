library('ROCR')
library('pROC')

if(!exists('hiv'))source('readData.R')

selector<-!apply(is.na(hiv[,goodTargetCols]),1,any)

pdf(file.path('out','roc.pdf'),width=4.5,height=4.5)
  for(genital in unique(hiv$isGenital)){
    message(ifelse(genital,'Genital','Plasma'))
    thisSelect<-selector&hiv$Selection=='UT'&(hiv$isGenital==genital|!hiv$isDonor)
    #calculate prediction and confidence interval
    preds<-lapply(goodTargetCols,function(xx)prediction(hiv[thisSelect,xx], !hiv[thisSelect,'isDonor']))
    sensCI<-lapply(goodTargetCols,function(xx)ci.se(roc(!hiv[thisSelect,'isDonor'],hiv[thisSelect,xx]),specificities=seq(0,1,.01),progress='text'))
    names(preds)<-goodTargetCols
    perfs<-lapply(preds,function(pred)performance(pred, "tpr", "fpr"))
    aucs<-lapply(preds,function(pred)performance(pred, "auc")@y.values[[1]])
    cols<-rainbow(length(perfs),alpha=.7)
    cols2<-rainbow(length(perfs),alpha=.3)
    names(cols)<-names(cols2)<-names(perfs)
    par(mar=c(3.5,3.5,1,.1))
    plot(1,1,type='n',xlim=0:1,ylim=0:1,xlab=perfs[[1]]@x.name,ylab=perfs[[1]]@y.name,las=1,mgp=c(2.2,.5,0),main=ifelse(genital,'Genital','Plasma'))
    #add ROC lines
    for(ii in 1:length(perfs))lines(perfs[[ii]]@x.values[[1]],perfs[[ii]]@y.values[[1]],col=cols[[ii]],lw=3)
    #add confidence intervals
    for(ii in 1:length(perfs))polygon(1-c(seq(0,1,.01),rev(seq(0,1,.01))),c(sensCI[[ii]][,'2.5%'],rev(sensCI[[ii]][,'97.5%'])),col=cols2[[ii]],border=NA)
    abline(0,1,lty=2)
    legend('bottomright',sub('\\(.*$','',names(perfs)),col=cols,lwd=2,adj=c(0,.5),bty='n')
  }
dev.off()
