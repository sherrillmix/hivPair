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



tf6<-readFaDir('gc','CH')
tf6$tf<-grepl('[tT]/?[fF]',tf6$name)
tf6<-tf6[order(tf6$file,!tf6$tf),]
tf6$seq<-toupper(tf6$seq)
write.fa(tf6$name,tf6$seq,'out/allTf6month.fa')
tf6Align<-readFaDir('gc/allAlign','fasta$')
for(ii in unique(tf6Align$file[!tf6Align$file %in% c('LTR.nt.fasta','AlignSeq.nt.fasta')])){
  seqs<-degap(substring(tf6Align[tf6Align$file==ii,'seq'],1,30))
  names(seqs)<-tf6Align[tf6Align$file==ii,'name']
  tf6[,sub('.nt.fasta','',ii)]<- mapply(function(xx,yy)gregexpr(xx,yy)[[1]],seqs[tf6$name],tf6$seq)
}
if(any(tf6[,sub('.nt.fasta','',unique(tf6Align$file[!tf6Align$file %in% c('LTR.nt.fasta','AlignSeq.nt.fasta')]))])==-1)stop('Problem finding proteins')


if(!file.exists('out/tf.png')){
  png('out/tf.png',height=3000,width=3000,res=300)
  plotDNA(tf6$seq)
  dev.off()
}
changes<-tapply(strsplit(tf6$seq,''),tf6$file,function(splits){
  diffs<-which(splits[[1]]!=splits[[2]])
  oneDown<-do.call(rbind,lapply(diffs,function(xx)substring(sapply(splits,paste,collapse=''),xx-1,xx)))
  oneUp<-do.call(rbind,lapply(diffs,function(xx)substring(sapply(splits,paste,collapse=''),xx,xx+1)))
  return(list(table(oneDown[,1],oneDown[,2]),table(oneUp[,1],oneUp[,2]),table(c(oneDown[,1],oneUp[,1]),c(oneDown[,2],oneUp[,2]))))
})

cgs<-sapply(names(changes),function(xx){
  message('# ',sub(' .*$','',xx),' #')
  up<-changes[[xx]][[1]]
  down<-changes[[xx]][[2]]
  message(' ',sum(up),' changes')
  message()
  zz<-changes[[xx]][[3]]
  twoMers<-as.vector(outer(c('A','C','T','G'),c('A','C','T','G'),paste,sep=''))
  changes<-do.call(c,lapply(twoMers,function(twoMer){
    out<-c(
      ifelse(twoMer %in% rownames(zz),sum(zz[twoMer,]),0),
      ifelse(twoMer %in% colnames(zz),sum(zz[,twoMer]),0)
    )
    names(out)<-sprintf(c('from%s','to%s'),twoMer)
    out
  }))
  return(changes)
})
colnames(cgs)<-sub(' .*$','',colnames(cgs))

determinants<-read.csv("gc/determinants.csv",stringsAsFactors=FALSE)
dets<-lapply(strsplit(determinants$Determinants.of.IFN.resistance, ', '),function(xx){
  xx<-xx[grepl('[0-9]',xx)]
  if(length(xx)==0)return(NULL)
  out<-do.call(rbind,strsplit(gsub('[()]','',xx),' '))
  if(any(probs<-grepl('-',out[,2]))){
    coords<-lapply(strsplit(sub('Δ','',out[probs,2]),'-'),as.numeric)
    reps<-lapply(coords,function(xx)xx[1]:xx[2])
    out<-rbind(out[!probs,], matrix(c(rep(out[probs,1],sapply(reps,length)),unlist(reps)),ncol=2))
  }
  return(out)
})

