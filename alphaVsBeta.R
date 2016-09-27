if(!exists('hiv'))source("readData.R")

cols<-rainbow.lab(length(unique(hiv$Pair.ID..)),alpha=.7)
names(cols)<-sort(unique(hiv$Pair.ID..))
pdf('out/alphaVsBeta.pdf')
  plot(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.,hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.,log='xy',xaxt='n',yaxt='n',xlab='IFNa2 IC50 (pg/ml)',ylab='IFNbeta IC50 (pg/ml)',bg=cols[as.character(hiv$Pair.ID..)],pch=21,col=NA)
  segments(min(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.),min(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.),max(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.),max(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.))
  logAxis(1)
  logAxis(2,las=1)
  legend('topleft',names(cols),pch=21,pt.bg=cols,inset=.01,title='Pair')
  for(ii in unique(hiv$select)){
    withAs(xx=hiv[hiv$select==ii,],plot(xx$IFNa2.Pooled.Donor.cells.IC50..pg..ml.,xx$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.,log='xy',xaxt='n',yaxt='n',xlab='IFNa2 IC50 (pg/ml)',ylab='IFNbeta IC50 (pg/ml)',bg=cols[as.character(xx$Pair.ID..)],pch=21,col=NA,main=ii)
    logAxis(1)
    logAxis(2,las=1)
  }
dev.off()

line<-function(alpha,beta){
  coefs<-lm(range(log10(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.))~range(log10(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.)))$coefficients
  pred<-coefs['(Intercept)']+coefs['range(log10(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.))']*log10(beta)
  return(log10(alpha)>pred)
}

scrambleWithinPatient<-function(alpha,pair){
  return(ave(alpha,pair,FUN=sample))
}

aboveLine<-mean(line(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.,hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.))
if(!exists('scrambleAbove'))scrambleAbove<-replicate(100000,mean(line(scrambleWithinPatient(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml.,hiv$Pair.ID..),hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml.)))
pdf('out/aboveLine.pdf')
  hist(scrambleAbove,xlim=c(.7,1))#range(c(aboveLine,scrambleAbove)))
  abline(v=aboveLine,col='red')
  axis(1,seq(.7,.95,.05))
dev.off()

