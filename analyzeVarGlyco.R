if(!exists('hiv'))source('readData.R')

donorMedians<-sapply(donorLengths,median)
lessMedianP<-sapply(donorLengths,function(x)mean(x<median(x)))
recLess<-unlist(recLengths)<rep(donorMedians,sapply(recLengths,length))
lessPs<-rep(lessMedianP,sapply(recLengths,length))

possibleObs<-do.call(expand.grid,rep(list(c(TRUE,FALSE)),length(lessPs)))
nLess<-apply(possibleObs,1,sum)
obsP<-apply(possibleObs,1,function(x)prod(lessPs[x])*prod((1-lessPs[!x])))
print(sum(obsP[nLess>=sum(recLess)]))

donorMedians<-sapply(donorGlyco,median)
lessMedianP<-sapply(donorGlyco,function(x)mean(x<median(x)))
recLess<-unlist(recGlyco)<rep(donorMedians,sapply(recGlyco,length))
lessPs<-rep(lessMedianP,sapply(recGlyco,length))

possibleObs<-do.call(expand.grid,rep(list(c(TRUE,FALSE)),length(lessPs)))
nLess<-apply(possibleObs,1,sum)
obsP<-apply(possibleObs,1,function(x)prod(lessPs[x])*prod((1-lessPs[!x])))
print(sum(obsP[nLess>=sum(recLess)]))


