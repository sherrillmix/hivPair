var<-read.csv('data/varLoops.csv',header=TRUE,check.names=FALSE)
colnames(var)<-sprintf('CH%s',colnames(var))
if(any(!colnames(var) %in% c(unlist(recs),unlist(donors))))stop(simpleError('Unknown donor/recipient'))
if(any(!c(unlist(recs),unlist(donors)) %in% colnames(var)))stop(simpleError('Missing donor/recipient'))

donorLengths<-lapply(donors,function(x)var[!is.na(var[,x]),x])
recLengths<-lapply(recs,function(x)lapply(x,function(y)var[!is.na(var[,y]),y]))
recLength<-lapply(recLengths,function(xx)sapply(xx,function(yy){
  lengthTab<-table(yy)
  maxLength<-names(lengthTab)[which.max(lengthTab)]
  if(any(lengthTab[maxLength]/sum(lengthTab)<.75))stop(simpleError('Unclear recipient length'))
  return(as.numeric(maxLength))
}))

glyco<-read.csv('data/glyco.csv',header=TRUE,check.names=FALSE)
colnames(glyco)<-sprintf('CH%s',colnames(glyco))
if(any(!colnames(glyco) %in% c(unlist(recs),unlist(donors))))stop(simpleError('Unknown donor/recipient'))
if(any(!c(unlist(recs),unlist(donors)) %in% colnames(glyco)))stop(simpleError('Missing donor/recipient'))

donorGlyco<-lapply(donors,function(x)glyco[!is.na(glyco[,x]),x])
recGlyco<-lapply(recs,function(x)lapply(x,function(y)glyco[!is.na(glyco[,y]),y]))
recGlyco<-lapply(recGlyco,function(xx)sapply(xx,function(yy){
  lengthTab<-table(yy)
  maxLength<-names(lengthTab)[which.max(lengthTab)]
  if(any(lengthTab[maxLength]/sum(lengthTab)<.5))stop(simpleError('Unclear recipient length'))
  return(as.numeric(maxLength))
}))


readLoopFile<-function(file){
  lines<-readLines(file)
  lines<-lines[lines!='']
  lines<-lines[!grepl('^Selected Region:',lines)]
  heads<-grep('Region\tName\t',lines)
  if(any(lines[heads]!=lines[heads][1]))stop(simpleError('Disagreement in headers'))
  lines<-lines[-heads[-1]]
  df<-read.table(textConnection(lines),header=TRUE,stringsAsFactors=FALSE)
  return(df)
}
condenseLoops<-function(lengths,regions,names){
  check<-tapply(regions,names,function(x){
    if(length(unique(x))!=length(x))stop(simpleError('Non unique regions'))
    sort(unique(x))
  })
  if(any(!sapply(check,function(x)all(x==check[[1]]))))stop(simpleError('Mixed up regions'))
  return(tapply(lengths,names,sum))
}
test<-readLoopFile('data/var_reg_char_16-08-18_6-51.txt')
test2<-condenseLoops(test$Length,test$Region,test$Name)
write.csv(data.frame('name'=names(test2),'length'=test2),'out/loopTest.csv',row.names=FALSE)
summary(sort(donorLengths[[1]])==sort(test2[grep('CH492',names(test2))]))
summary(sort(recLengths[[1]][[1]])==sort(test2[grep('CH427',names(test2))]))


phyloDist<-data.frame('dist'=t(read.csv('data/median TF dists.csv')))
bDist<-read.csv('data/combinedB.3prUTpl.161020exc.dists.csv',header=FALSE,stringsAsFactors=FALSE)[,-1:-2]
cDist<-read.csv('data/combinedC.3prUTpl.161020exc.dists.csv',header=FALSE,stringsAsFactors=FALSE)[,-1:-2]
colnames(bDist)<-colnames(cDist)<-c('compare','name','dist')
allDist<-rbind(cDist,bDist)
allDist<-allDist[allDist$compare %in% c('CONSENSUS_B','CONSENSUS_C'),]
