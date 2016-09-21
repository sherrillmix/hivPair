library('rstan')

#for parallel
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores()

if(!exists('hiv'))source('readData.R')

targetCol<-'IFNa2.Pooled.Donor.cells.IC50..pg..ml.'
stanCode<-"
  data {
    int<lower=0> N;
    real ic50[N];
    int<lower=0,upper=1> isCladeB[N];
    int<lower=0,upper=1> isGenital[N];
    int<lower=0,upper=1> isRecipient[N];
    int<lower=0> nGroup;
    int<lower=1> group[N];
    int<lower=0> nPair;
    int<lower=1> pair[N];
  }
  parameters {
    real mus[nPair];
    real<lower=0.0000001> sigmas[nGroup];
    real genital;
    real recipient[nPair];
    real clade;
  }
  model {
    for (ii in 1:N)ic50[ii]~normal(
        mus[pair[ii]]+genital*isGenital[ii]+recipient[pair[ii]]*isRecipient[ii]+clade*isCladeB[ii]
        ,sigmas[group[ii]]
      );
  }
"

dat<-list(ic50=log10(hiv[,targetCol]),N=nrow(hiv),isCladeB=as.integer(hiv$Subtype=='B'),isGenital=as.integer(hiv$fluid=='PL'),isRecipient=as.integer(!hiv$donor),nGroup=max(hiv$group),group=hiv$group,pair=hiv$Pair.ID..,nPair=max(hiv$Pair.ID..))
fit <- stan(model_code = stanCode, data = dat, iter=1000, chains=4)


