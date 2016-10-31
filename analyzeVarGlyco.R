if(!exists('hiv'))source('readData.R')

donorMedians<-sapply(donorLengths,median)
lessMedianP<-sapply(donorLengths,function(x)mean(x<median(x)))
recLess<-unlist(recLength)<rep(donorMedians,sapply(recLength,length))
lessPs<-rep(lessMedianP,sapply(recLength,length))

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


quantP<-function(x,singleY)mean(singleY>=c(x,singleY))
fisherMethod<-function(ps)1-pchisq(-2*sum(log(ps)),length(ps)*2)

tfs<-phyloDist[grep("TF",rownames(phyloDist)),,drop=FALSE]
tfs$base<-sub('(v[0-9])?\\.TF','',rownames(tfs))
tfs$donor<-sapply(tfs$base,function(xx)hiv[hiv$donor & hiv$Pair.ID==unique(hiv[hiv$baseName==xx,'Pair.ID']),'baseName'][1])
allDist$base<-sapply(strsplit(allDist$name,'_'),'[[',1)
tfs$ps<-sapply(rownames(tfs),function(xx)quantP(allDist[allDist$base==tfs[xx,'donor'],'dist'],tfs[xx,'dist']))
fisherMethod(tfs$ps)



