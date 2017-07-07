if(!exists('rawSeqMat'))source('parseSeqs.R')

allNonGap<-range(which(apply(rawSeqMat!='-',2,all)))
alignedSeqs<-apply(rawSeqMat[,allNonGap[1]:allNonGap[2]],1,paste,collapse='')
degappedSeqs<-degap(alignedSeqs)
gcPos<-gregexpr('CG',degappedSeqs)
gcs<-sapply(gcPos,length)/nchar(degappedSeqs)
gcaPos<-gregexpr('[AT]CG[AT]',degappedSeqs)
gcas<-sapply(gcaPos,length)/nchar(degappedSeqs)
gcaPerGc<-gcas/gcs
pdf('out/gc.pdf')
  par(mgp=c(3,.2,0),tcl=-.1)
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcs[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='CG/nucleotide',las=1)
  for(ii in unique(hiv$Pair.ID)){
    main<-paste(unique(hiv$baseName[hiv$Pair.ID==ii][order(!hiv$donor[hiv$Pair.ID==ii])]),collapse=',')
    patCh<-20+as.numeric(as.factor(hiv$baseName[hiv$Pair.ID==ii&hiv$select=='UT']))
    vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$Pair.ID==ii&hiv$select=='UT'],gcs[hiv$Pair.ID==ii&hiv$select=='UT'],pch=patCh,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID][hiv$Pair.ID==ii],ylab='CG/nucleotide',main=main,las=1,ylim=range(gcs))
  }
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcas[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='[AT]CG[AT]/nucleotide',las=1)
  for(ii in unique(hiv$Pair.ID)){
    main<-paste(unique(hiv$baseName[hiv$Pair.ID==ii][order(!hiv$donor[hiv$Pair.ID==ii])]),collapse=',')
    patCh<-20+as.numeric(as.factor(hiv$baseName[hiv$Pair.ID==ii&hiv$select=='UT']))
    vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$Pair.ID==ii&hiv$select=='UT'],gcas[hiv$Pair.ID==ii&hiv$select=='UT'],pch=patCh,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID][hiv$Pair.ID==ii],ylab='[AT]CG[AT]/nucleotide',main=main,las=1,ylim=range(gcas))
  }
  vpPlot(ifelse(hiv$donor,'Donor','Recipient')[hiv$select=='UT'],gcaPerGc[hiv$select=='UT'],pch=21,bg=rainbow(max(hiv$Pair.ID))[hiv$Pair.ID[hiv$select=='UT']],ylab='[AT]CG[AT]/CG',las=1)
  plotDNA(noGapSeqs)
dev.off()


refCoords<-refPos[allNonGap[1]:allNonGap[2]]
prettyPos<-pretty(refCoords)
prettyPos[prettyPos>max(refCoords)]<-max(refCoords)
prettyPos[prettyPos<min(refCoords)]<-min(refCoords)
prettyX<-sapply(prettyPos,function(xx)min(which(refCoords==xx)))
gcPos<-gregexpr('C-*G',alignedSeqs)
gcaPos<-gregexpr('[AT]-*C-*G-*[AT]',alignedSeqs)
png('out/gcPos.png',height=4000,width=3000,res=250)
  par(mar=c(4,4,.2,4.5))
  plot(1,1,type='n',xlim=range(unlist(gcPos)),ylim=c(1-1,length(gcPos)+1),ylab='Virus',yaxs='i',las=1,xlab='Position (HXB2)',xaxt='n')
  axis(1,prettyX,prettyPos)
  abline(h=1:length(gcPos),col='#00000033')
  ordering<-rank(paste(hiv$Pair.ID,hiv$donor,hiv$Renamed))
  points(unlist(gcPos),rep(ordering,sapply(gcPos,length)),cex=.4)
  points(unlist(gcaPos),rep(ordering,sapply(gcaPos,length)),cex=.4,col='red')
  bottom<-sort(tapply(ordering,hiv$baseName,FUN=min))
  top<-tapply(ordering,hiv$baseName,FUN=max)[names(bottom)]
  labPos<-par('usr')[2]+diff(par('usr')[1:2])*rep(c(.01,.015),length.out=length(top))
  segments(labPos,bottom,labPos,top,xpd=NA)
  text(labPos+diff(par('usr')[1:2])*.01,(top+bottom)/2,names(bottom),xpd=NA,adj=0)
dev.off()
png('out/gcPos2.png',height=4000,width=3000,res=250)
  par(mar=c(4,4,.2,4.5))
  plot(1,1,type='n',xlim=range(unlist(gcPos)),ylim=c(1-1,length(gcPos)+1),ylab='Virus',yaxs='i',las=1,xlab='Position (HXB2)',xaxt='n')
  axis(1,prettyX,prettyPos)
  abline(h=1:length(gcPos),col='#00000033')
  ordering<-rank(paste(hiv$donor,hiv$Pair.ID,hiv$Renamed))
  points(unlist(gcPos),rep(ordering,sapply(gcPos,length)),cex=.4)
  points(unlist(gcaPos),rep(ordering,sapply(gcaPos,length)),cex=.4,col='red')
  bottom<-sort(tapply(ordering,hiv$baseName,FUN=min))
  top<-tapply(ordering,hiv$baseName,FUN=max)[names(bottom)]
  labPos<-par('usr')[2]+diff(par('usr')[1:2])*rep(c(.01,.015),length.out=length(top))
  segments(labPos,bottom,labPos,top,xpd=NA)
  text(labPos+diff(par('usr')[1:2])*.01,(top+bottom)/2,names(bottom),xpd=NA,adj=0)
dev.off()

for(ii in unique(hiv$Pair.ID)){
  pats<-unique(hiv$baseName[hiv$Pair.ID==ii][order(!hiv$donor[hiv$Pair.ID==ii])])
  png(sprintf('out/gcPos_%s.png',paste(pats,collapse='-')),height=2000,width=2500,res=250)
    par(mar=c(4,5,.1,.1))
    selector<-hiv$Pair.ID==ii&hiv$select=='UT'&hiv$fluid=='PL'
    plot(1,1,type='n',xlim=range(unlist(gcPos))+c(-20,20),ylim=c(1-1,length(gcPos[selector])+1),ylab='',yaxs='i',las=1,xlab='Position (HXB2)',yaxt='n',xaxt='n',xaxs='i')
    axis(1,prettyX,prettyPos)
    abline(h=1:length(gcPos),col=ifelse(hiv[selector,'donor'][order(hiv[selector,'donor'])],'#0000FF99','#00000033'))
    points(unlist(gcPos[selector][order(hiv[selector,'donor'])]),rep((1:length(gcPos[selector])),sapply(gcPos[selector],length)),cex=.6)
    points(unlist(gcaPos[selector][order(hiv[selector,'donor'])]),rep((1:length(gcaPos[selector])),sapply(gcaPos[selector],length)),cex=.6,col='red')
    axis(2,1:sum(selector),hiv[selector,'Renamed'][order(hiv[selector,'donor'])],cex.axis=.5,las=1,mgp=c(1,.15,0),tcl=-.05)
  dev.off()
}
