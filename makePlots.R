if(!exists('hiv'))source('readData.R')



varCols<-which(colnames(hiv)=='Replicative.capacity.Single.Donor.p24.d7'):which(colnames(hiv)=='Bnaber.IC50')
pdf('out/pairs.pdf',width=12,height=12)
thisData<-hiv[,varCols]
colnames(thisData)<-gsub('\\.+','\n',colnames(thisData))
plot(thisData,col='#00000077',pch='.',cex=2.5)
dev.off()

#align<-levenAlign(hiv$seq,hiv$seq[1])



uniqSample<-unique(hiv$sample[order(hiv$Pair.ID..,hiv$Donor.or.Recipient)])
pairId<-sub('.* ','',uniqSample)
labs<-factor(hiv$sample,levels=uniqSample)
cols<-c('#FF000066','#0000FF66')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
changes<-which(isChange)
labPos<-1:length(uniqSample)+cumsum(isChange)*.5
xPos<-labPos[as.numeric(labs)]
pdf('out/untreated_vpPlot.pdf',height=4,width=8)
for(ii in varCols){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)&hiv$select=='UT'&hiv$fluid=='PL'
  par(mar=c(6.5,4,.1,.1),mgp=c(3,.7,0))
  plot(1,1,las=2,ylab=colnames(hiv)[ii],xlim=range(labPos),ylim=range(thisData,na.rm=TRUE),xaxt='n',xlab='',type='n')
  axis(1,labPos,uniqSample,las=2)
  xOffset<-offsetX(thisData[selector],labs[selector])
  rect(c(1,labPos[changes])-.75,par('usr')[3],c(labPos[changes],length(pairId)+1)-.75,par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=.9,col=NA,bg=cols[hiv$donor[selector]+1],pch=21)
}
dev.off()






uniqSample<-unique(hiv$sampleSelect[order(hiv$Pair.ID..,hiv$donorRec,c('UT'=1,'A2'=2,'BE'=3)[hiv$select])])
pairId<-sub('^[A-Za-z0-9-]+ ([0-9]) [A-Z0-9]+','\\1',uniqSample)
labs<-factor(hiv$sampleSelect,levels=uniqSample)
cols<-c('#FF000066','#0000FF66')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
donorRec<-sub(' .+','',uniqSample)
treat<-sub('.* ','',uniqSample)
isChange2<-c(FALSE,donorRec[-1]!=donorRec[-length(pairId)])
changes<-which(isChange)
changes2<-which(c(FALSE,treat[-1]!=treat[-length(pairId)]))
labPos<-1:length(uniqSample)+cumsum(isChange)*.5+cumsum(isChange2*.5)
xPos<-labPos[as.numeric(labs)]
cols<-rainbow.lab(10,alpha=.7)[c(1:3,8:10)]
names(cols)<-unique(paste(hiv$donor,hiv$select))

pdf('out/selects_vpPlot.pdf',height=4,width=8)
for(ii in varCols){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)
  par(mar=c(7.3,4,.3,.1),mgp=c(3,.7,0))
  plot(1,1,las=2,ylab=colnames(hiv)[ii],xlim=range(labPos)+c(-.7,.7),ylim=range(thisData,na.rm=TRUE),xaxt='n',xlab='',type='n',log='y',xaxs='i')
  axis(1,labPos,uniqSample,las=2)
  xOffset<-offsetX(log10(thisData[selector]),labs[selector],width=.4)
  #xOffset<-ave(thisData[selector],labs[selector],FUN=function(x)swarmx(rep(0,length(x)),x,log='y',cex=.9)$x)
  rect(c(1,labPos[changes])-1,10^par('usr')[3],c(labPos[changes],length(pairId)+1)-1,10^par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  #rect(c(1,labPos[changes2])-.5,10^par('usr')[3],c(labPos[changes2],length(pairId)+1)-.5,10^par('usr')[4],col=rep(c(NA,'#00000006'),length.out=length(changes2)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=.9,col=NA,bg=cols[paste(hiv$donor,hiv$select)[selector]],pch=21)
}
dev.off()


#Env  content ( Env :RT ratios)
#Infectivity
#Replicative capacity
#IFN alpha IC50
selectVars<-c(
  'Env.RT'='Env/RT',
  'Infectivity.RLU.pg.RT...T1249.'='Infectivity (RLU/pg RT)',
  'Replicative.capacity.Pool.Donor.p24.d7'='Pooled donor\nReplicative capactity (day 7 p24)',
  'IFNa2.PD.IC50..U.ml.'='Pooled donor\nIFNa2 IC50 (U/ml)'
)
  #'Replicative.capacity.Single.Donor.p24.d7'='Single donor\nReplicative capactity\n(day 7 p24)',
  #'IFNa2.SD.IC50..U.ml.'='Single donor\nIFNa2 IC50 (U/ml)',

logAxis<-function(x,axisNum=2,addExtra=TRUE,spreadRange=1.3,...){
  minX<-min(log10(x),na.rm=TRUE) 
  maxX<-max(log10(x),na.rm=TRUE) 
  allTicks<-unlist(lapply(floor(minX):ceiling(maxX),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<10^maxX*spreadRange & allTicks>10^minX/spreadRange]
  axis(axisNum,allTicks,rep('',length(allTicks)),tcl=-.2)
  prettyY<-seq(ceiling(log10(min(x,na.rm=TRUE))),floor(log10(max(x,na.rm=TRUE))),1)
  if(addExtra){
    if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(5),prettyY-log10(2)))
    if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(2),prettyY-log10(5)))
  }
  axis(axisNum,10^prettyY,10^prettyY,las=2,...)
}
uniqSample<-unique(hiv$sampleFluid[order(hiv$Pair.ID..,hiv$Donor.or.Recipient,c('PL'=1,'CV'=2,'SE'=3)[hiv$fluid])])
pairId<-sub('.* ([0-9]+) [A-Z]+','\\1',uniqSample)
labs<-factor(hiv$sampleFluid,levels=uniqSample)
cols<-c('TRUE PL'='#FF000055','TRUE CV'='#0000FF55','TRUE SE'='#0000FF55','FALSE PL'='#FF770055')
isChange<-c(FALSE,pairId[-1]!=pairId[-length(pairId)])
changes<-which(isChange)
labPos<-1:length(uniqSample)+cumsum(isChange)*.5
xPos<-labPos[as.numeric(labs)]
pdf('out/selectVar.pdf',height=8,width=7)
layout(matrix(c(length(selectVars)+2,1:(length(selectVars)+1)),ncol=1),height=c(.05,rep(1,length(selectVars)),.6))
for(ii in names(selectVars)){
  thisData<-hiv[,ii]
  selector<-!is.na(thisData)&hiv$select=='UT'
  par(mar=c(0,4.8,0,.2),mgp=c(3,.7,0),lheight=.8)
  plot(1,1,las=2,ylab=selectVars[ii],xlim=range(labPos),ylim=range(thisData,na.rm=TRUE)*c(.8,1.2),xaxt='n',xlab='',type='n',log='y',yaxt='n')
  logAxis(thisData)
  xOffset<-offsetX(thisData[selector],labs[selector],width=.45,varwidth=TRUE)
  rect(c(1,labPos[changes])-.75,10^par('usr')[3],c(labPos[changes],length(pairId)+1)-.75,10^par('usr')[4],col=rep(c(NA,'#00000011'),length.out=length(changes)),border=NA)
  points(xPos[selector]+xOffset,thisData[selector],cex=1,col=NA,bg=cols[paste(hiv$donor,hiv$fluid)[selector]],pch=21)
  #axis(1,labPos,rep('',length(labPos)))
}
axis(1,labPos,uniqSample,las=2)
dev.off()






for(fluid in list('PL',unique(hiv$fluid))){
  for(plotType in c('box','points')){
    for(var in names(selectVars)){
      message(var)
      thisData<-hiv[hiv$select=='UT'&hiv$fluid %in% fluid&!is.na(hiv[,var]),]
      infectData<-tapply(thisData[,var][thisData$donor],thisData$Pair.ID..[thisData$donor],median,na.rm=TRUE)
      #infectData<-sort(infectData,decreasing=TRUE)
      ylim<-range(thisData[,var],na.rm=TRUE)
      pairOrder<-order(infectData)
      pairOrder<-1:length(infectData)
      cols<-rainbow.lab(length(infectData),alpha=.5)
      cols2<-rainbow.lab(length(infectData))
      names(cols)<-names(cols2)<-sort(unique(thisData$Pair.ID..))
      #rectWidth<-.004
      nDonor<-length(infectData)
      donorStep<-.6/(nDonor-1)
      nRecs<-length(unique(paste(thisData$Donor.or.Recipient,thisData$Pair.ID..)[!thisData$donor]))
      recStep<-.6/(nRecs-1)
      pdf(sprintf('out/pair/pair_%s_%s_%s.pdf',plotType,var,paste(fluid,collapse='-')),width=4,height=4)
        par(mar=c(6.1,3,.1,.1))
        plot(1,1,type='n',xlim=c(.6,2.4),ylim=ylim,ylab=selectVars[var],xlab='',xaxt='n',log='y',las=1,mgp=c(2,.7,0),yaxt='n')
        logAxis(ylim,mgp=c(3,.8,0),addExtra=FALSE)
        axis(1,1:2,c('Donor','Recipient'),mgp=c(3,.1,0),tcl=0)
        box()
        #abline(v=1.5)
        donorPos<--.3
        recPos<--.3
        pointSize<-.4
        for(ii in 1:length(infectData)){
          thisPair<-names(infectData)[ii]
          thisDonor<-thisData[thisData$Pair.ID..==thisPair&thisData$donor,]
          thisRec<-thisData[thisData$Pair.ID..==thisPair&!thisData$donor,]
          thisRec<-split(thisRec,thisRec$Donor.or.Recipient)
          dPos<-1+donorPos
          rPos<-2+recPos+0:(length(thisRec)-1)*recStep
          for(jj in 1:length(thisRec)){
            #thisP<-wilcox.test(thisRec[[jj]][,var],thisDonor[,var])$p.value
            #thisCol<-ifelse(thisP<.05,cols[thisPair],gray(0,alpha=.1))
            thisCol<-cols[thisPair]#ifelse(thisP<.05,cols[thisPair],gray(0,alpha=.1))
            if(plotType=='points'){
              segments(dPos,median(thisDonor[,var]),rPos[[jj]],median(thisRec[[jj]][,var]),col=thisCol,lty=1,lwd=2)
              segments(rPos[[jj]],min(thisRec[[jj]][,var]),rPos[[jj]],max(thisRec[[jj]][,var]),col=cols2[thisPair])
              offsetPos<-offsetX(thisRec[[jj]][,var],width=.025)
              #offsetPos<-swarmx(rep(0,nrow(thisRec[[jj]])),thisRec[[jj]][,var],cex=.4,log='y')$x
              points(rep(rPos[[jj]],nrow(thisRec[[jj]]))+offsetPos,thisRec[[jj]][,var],bg=cols[thisPair],pch=21,cex=pointSize,col=NA)
              segments(rPos[[jj]]+.02,median(thisRec[[jj]][,var]),rPos[[jj]]-.02,median(thisRec[[jj]][,var]),col=cols[thisPair])
            }else if(plotType=='box'){
              segments(dPos,median(thisDonor[,var]),rPos[[jj]],median(thisRec[[jj]][,var]),col=thisCol,lty=1,lwd=2)
              segments(rPos[[jj]],min(thisRec[[jj]][,var]),rPos[[jj]],max(thisRec[[jj]][,var]),col=cols2[thisPair])
              box<-boxplot(thisRec[[jj]][,var],plot=FALSE)
              #rect(rPos[[jj]]+.02,box$conf[1,1],rPos[[jj]]-.02,box$conf[2,1],col=cols2[thisPair])
              rect(rPos[[jj]]+.02,box$stats[2,1],rPos[[jj]]-.02,box$stats[4,1],col=cols2[thisPair])
              segments(rPos[[jj]]+.02,median(thisRec[[jj]][,var]),rPos[[jj]]-.02,median(thisRec[[jj]][,var]))
            }
          }
          if(plotType=='points'){
            segments(dPos,min(thisDonor[,var]),dPos,max(thisDonor[,var]),col=cols2[thisPair])
            offsetPos<-offsetX(thisDonor[,var],width=.025)
            #offsetPos<-swarmx(rep(0,nrow(thisDonor)),thisDonor[,var],cex=.4,log='y')$x
            points(rep(dPos,nrow(thisDonor))+offsetPos,thisDonor[,var],bg=cols[thisPair],pch=21,cex=pointSize,col=NA)
            segments(dPos+.02,median(thisDonor[,var]),dPos-.02,median(thisDonor[,var]),col=cols[thisPair])
          }else if(plotType=='box'){
            box<-boxplot(thisDonor[,var],plot=FALSE)
            segments(dPos,min(thisDonor[,var]),dPos,max(thisDonor[,var]),col=cols2[thisPair])
            #rect(dPos+.02,box$conf[1,1],dPos-.02,box$conf[2,1],col=cols2[thisPair])
            rect(dPos+.02,box$stats[2,1],dPos-.02,box$stats[4,1],col=cols2[thisPair])
            segments(dPos+.02,median(thisDonor[,var]),dPos-.02,median(thisDonor[,var]))
          }
          donorPos<-donorPos+donorStep
          recPos<-recPos+recStep*length(thisRec)
        }
        #legend('bottomright',sprintf('%s ',names(cols)),lwd=2,col=cols,pt.cex=.4,pt.bg=cols2,ncol=2,x.intersp=.2,inset=.02,title='Pair',bty='o',cex=.9,box.col='#00000055') #,pch=21
        legend(1.5,convertLineToUser(1.5,1),sprintf('%s',pairNames[names(cols)]),lwd=2,col=cols2,pt.cex=.4,pt.bg=cols2,ncol=2,inset=-.02,title='Pair',bty='o',cex=.7,xpd=NA,yjust=1,xjust=.5) #,pch=21
        #legend('bottomright',sprintf('%s ',names(cols)),lwd=2,col=cols,pt.cex=.4,pt.bg=cols2,ncol=2,x.intersp=.2,inset=.02,title='Pair',bty='o',cex=.9,box.col='#00000055') #,pch=21
      dev.off()
      pdf(sprintf('out/pair/box_%s.pdf',var),width=7,height=4)
       par(las=2,mar=c(6,4,.1,.1))
       boxplot(thisData[,var]~paste(thisData$Pair.ID..,thisData$Donor.or.Recipient),ylab=selectVars[var],log='y',col=unlist(lapply(nRecs,function(x)rep(c('red','blue'),c(1,x)))))
      dev.off()
    }
  }
}

