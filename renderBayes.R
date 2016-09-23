
if(!exists('hiv'))source('readData.R')
tmp<-tempfile()

out<-'
digraph G {
  node [shape=plaintext fontname="Arial"];
  {rank=source donor genital clade recipient}
  {rank=sink varGenital varDonor varRecipient}
  {rank=same _DISTS_HERE_}

  donor [label="Donor",shape="ellipse"];
  genital [label="Genital",shape="ellipse"];
  clade [label="CladeB",shape="ellipse"];
  recipient [label="Recipient",shape="ellipse"];
'

dists<-c()
vars<-c()
for(ii in sort(unique(hiv$Pair.ID..))[c(1,2)]){
  thisData<-hiv[hiv$Pair.ID..==ii,]
  dists<-c(dists,sprintf(c('pair%dRecipientDist','pair%dDonorDist'),ii))
  vars<-c(vars,sprintf(c('varDonor_%d','varRecipient_%d'),ii))
  insert<-'
  recipient  -> recipient_ID_HERE_;
  donor -> pair_ID_HERE_donor;
  recipient_ID_HERE_ -> pair_ID_HERE_recipient [dir="none"];
  pair_ID_HERE_donor -> pair_ID_HERE_recipient [dir="none"];
  pair_ID_HERE_recipient -> pair_ID_HERE_RecipientDist;
  pair_ID_HERE_donor -> pair_ID_HERE_DonorDist;
  varDonor_ID_HERE_ -> pair_ID_HERE_DonorDist;
  varDonor -> varDonor_ID_HERE_;
  varRecipient -> varRecipient_ID_HERE_;
  varRecipient_ID_HERE_ -> pair_ID_HERE_RecipientDist;
  recipient_ID_HERE_ [label="Recipient_ID_HERE_ Effect",shape="rectangle"];
  pair_ID_HERE_donor [label="Pair _ID_HERE_ Donor Mean",shape="rectangle"];
  pair_ID_HERE_recipient [style=invis, height=0, label=""];
  varRecipient_ID_HERE_ [label="Recipient _ID_HERE_ Variance",shape="rectangle"];
  varDonor_ID_HERE_ [label="Donor _ID_HERE_ Variance",shape="rectangle"];
  pair_ID_HERE_DonorDist [label="Pair _ID_HERE_ Donor",shape="oval"];
  pair_ID_HERE_RecipientDist [label="Pair _ID_HERE_ Recipient",shape="oval"];
  '
  if(any(thisData$fluid!='PL')){
    dists<-c(dists,sprintf('pair%dGenitalDist',ii))
    vars<-c(vars,sprintf(c('varGenital_%d'),ii))
    insert<-c(insert,'
      genital  -> genital_ID_HERE_;
      genital_ID_HERE_ -> pair_ID_HERE_genital [dir="none"];
      pair_ID_HERE_donor -> pair_ID_HERE_genital [dir="none"];
      pair_ID_HERE_genital -> pair_ID_HERE_GenitalDist;
      varGenital_ID_HERE_ -> pair_ID_HERE_GenitalDist;
      varGenital -> varGenital_ID_HERE_;
      genital_ID_HERE_ [label="Genital_ID_HERE_ Effect",shape="rectangle"];
      varGenital_ID_HERE_ [label="Genital _ID_HERE_ Variance",shape="rectangle"];
      pair_ID_HERE_GenitalDist [label="Pair _ID_HERE_ Genital",shape="oval"];
      pair_ID_HERE_genital [style=invis, height=0, label=""];
    ')
    genIc<-sprintf('genitalIC50_%d_%d',ii,1:sum(thisData$fluid!='PL'))
    genIc50s<-sprintf('pair%dGenitalDist -> %s;\n%s [label="IC50",shape="diamond",fontsize=6,width=.3,fixedsize=TRUE];',ii,genIc,genIc);
    out<-c(out,genIc50s)
  }
  donorIc<-sprintf('donorIC50_%d_%d',ii,1:sum(thisData$fluid=='PL'&thisData$donor))
  recIc<-sprintf('recipientIC50_%d_%d',ii,1:sum(thisData$fluid=='PL'&!thisData$donor))
  recIc50s<-sprintf('pair%dRecipientDist-> %s;\n%s [label="IC50",shape="diamond",fontsize=6,width=.3,fixedsize=TRUE];',ii,recIc,recIc);
  donorIc50s<-sprintf('pair%dDonorDist -> %s;\n%s [label="IC50",shape="diamond",fontsize=6,width=.3,fixedsize=TRUE];',ii,donorIc,donorIc);
  out<-c(out,recIc50s,donorIc50s)
  if(thisData$Subtype[1]=='B'){
    insert<-c(insert,'
      clade  -> clade_ID_HERE_;
      clade_ID_HERE_ -> pair_ID_HERE_recipient [dir="none"];
      clade_ID_HERE_ [label="CladeB_ID_HERE_ Effect",shape="rectangle"];
    ')
  }
  insert<-gsub('_ID_HERE_',ii,insert)
  out<-c(out,insert)
}
  #pair1DonorDist -> ic50_1;
  #pair1DonorDist -> ic50_2;
  #pair1GenitalDist -> ic50_3;
  #ic50_1 [label="IC50",shape="diamond"];
  #ic50_2 [label="IC50",shape="diamond"];
  #ic50_3 [label="IC50",shape="diamond"];


out<-sub('_DISTS_HERE_',paste(dists,collapse=' '),out)
out<-c(out,'
  varRecipient [label="Recipient Variance",shape="ellipse"];
  varDonor [label="Donor Variance",shape="ellipse"];
  varGenital [label="Genital Variance",shape="ellipse"];
}')


writeLines(out,tmp)
system(sprintf('dot -Tpng %s >out/bayesModel.png',tmp))
