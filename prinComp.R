if(!exists('hiv'))source('readData.R')

targetCols<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pooled.Donor.cells.p24.d7'='Pooled donor\nReplicative capacity (day 7 p24)',
  'IFNbeta.Pooled.Donor.cells.IC50..pg.ml.'='IFNbeta IC50 (pg/ml)',
  'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'='IFNa2 IC50 (pg/ml)'
)

selector<-!apply(is.na(hiv[,names(targetCols)]),1,any)
pca<-prcomp(hiv[selector,names(targetCols)],scale.=TRUE)
pcaPoints<-t(t(pca$x[,1:2])/pca$sdev[1:2]/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
importance<-summary(pca)$importance[2,]

cols<-rainbow.lab(length(unique(hiv$donor)),alpha=.8)
names(cols)<-unique(hiv$donor)
pch<-structure(c(21,22,22),names=c('PL','SE','CV'))
col=ifelse(hiv[selector,'select']=='BE','black',ifelse(hiv[selector,'select']=='A2','red','#00000055'))

pdf('out/pca.pdf')
  plot(
    1,1,type='n',las=1,
    xlim=range(pcaPoints[,1]),ylim=range(pcaPoints[,2]),
    xlab=sprintf('Principal component 1 (%d%% of variance)',round(importance[1]*100)),
    ylab=sprintf('Principal component 2 (%d%% of variance)',round(importance[2]*100))
  )
  for(ii in c('UT','A2','BE')){
    selectSelector<-hiv[selector,'select']==ii
    points(pcaPoints[selectSelector,],bg=cols[as.character(hiv[selector,'donor'])][selectSelector],pch=pch[hiv[selector,'fluid']][selectSelector],col=col[selectSelector])
  }
  #plot(pca)
dev.off()
