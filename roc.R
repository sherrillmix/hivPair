library(ROCR)
library(pROC)

library(doParallel)
registerDoParallel(cl <- makeCluster(getOption("mc.cores", parallel::detectCores())))
pdf('out/roc.pdf',width=4.5,height=4.5)
  for(genital in unique(hiv$isGenital)){
    message(ifelse(genital,'Genital','Plasma'))
    thisSelect<-selector&hiv$select=='UT'&(hiv$isGenital==genital|!hiv$donor)
    preds<-lapply(names(goodTargetCols),function(xx)prediction(hiv[thisSelect,xx], !hiv[thisSelect,'donor']))
    aucCI<-lapply(names(goodTargetCols),function(xx)ci(!hiv[thisSelect,'donor'],hiv[thisSelect,xx]))
    sensCI<-lapply(names(goodTargetCols),function(xx)ci.se(roc(!hiv[thisSelect,'donor'],hiv[thisSelect,xx]),specificities=seq(0,1,.01),parallel=TRUE))
    names(preds)<-names(goodTargetCols)
    perfs<-lapply(preds,function(pred)performance(pred, "tpr", "fpr"))
    aucs<-lapply(preds,function(pred)performance(pred, "auc")@y.values[[1]])
    print(aucs)
    print(aucCI)
    cols<-rainbow.lab(length(perfs),alpha=.7)
    cols2<-rainbow.lab(length(perfs),alpha=.3)
    par(mar=c(3.5,3.5,1,.1))
    plot(1,1,type='n',xlim=0:1,ylim=0:1,xlab=perfs[[ii]]@x.name,ylab=perfs[[ii]]@y.name,las=1,mgp=c(2.2,.5,0),main=ifelse(genital,'Genital','Plasma'))
    for(ii in 1:length(perfs))lines(perfs[[ii]]@x.values[[1]],perfs[[ii]]@y.values[[1]],col=cols[[ii]],lw=3)
    for(ii in 1:length(perfs))polygon(1-c(seq(0,1,.01),rev(seq(0,1,.01))),c(sensCI[[ii]][,'2.5%'],rev(sensCI[[ii]][,'97.5%'])),col=cols2[[ii]],border=NA)
    abline(0,1,lty=2)
    legend('bottomright',sub('\\(.*$','',targetCols[names(perfs)]),col=cols,lwd=2,adj=c(0,.5),bty='n')
  }
  for(treat in unique(hiv$select[hiv$select!='UT'])){
    thisSelect<-selector&((hiv$donor&hiv$select==treat)|(!hiv$donor&hiv$select=='UT'))&!hiv$isGenital
    preds<-lapply(names(goodTargetCols),function(xx)prediction(hiv[thisSelect,xx], !hiv[thisSelect,'donor']))
    names(preds)<-names(goodTargetCols)
    perfs<-lapply(preds,function(pred)performance(pred, "tpr", "fpr"))
    aucs<-lapply(preds,function(pred)performance(pred, "auc")@y.values[[1]])
    cols<-rainbow.lab(length(perfs),alpha=.7)
    par(mar=c(3.5,3.5,1,.1))
    plot(1,1,type='n',xlim=0:1,ylim=0:1,xlab=perfs[[ii]]@x.name,ylab=perfs[[ii]]@y.name,las=1,mgp=c(2.2,.5,0),main=sprintf('%s donor vs UT recipient',treat))
    for(ii in 1:length(perfs))lines(perfs[[ii]]@x.values[[1]],perfs[[ii]]@y.values[[1]],col=cols[[ii]],lw=3)
    abline(0,1,lty=2)
    legend('bottomright',sub('\\(.*$','',targetCols[names(perfs)]),col=cols,lwd=2,adj=c(0,.5),bty='n')
  }
dev.off()
stopCluster(cl)


withAs(xx=hiv[hiv$select!='UT'&selector&!hiv$isGenital,],vpPlot(xx$donor,xx$IFNbeta.Pooled.Donor.cells.IC50..pg.ml))