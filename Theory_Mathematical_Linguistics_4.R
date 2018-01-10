library(protr)
library(CHNOSZ)
library(bio3d)
library(VarfromPDB)
library(dplyr)
library(HMM)
library(markovchain)
library(mcmc)
library(combinat)
library(permute)
library(permutations)
library(Peptides)
library(boot)
library(stringi)
library(stringr)
library(utils)
library(DescTools)
library(seqinr)
library(xtable)

#----------------------------------------------Data--------------------------------------------------------
data(thermo)
data(aa.table)
data(aa.index)
data(aa)
data(pepdata)
peptide.properties<-colnames(pepdata)

Amino.Study<-c("A","P","S","T")
amino.study.df<-data.frame()

Amino.Study.Notes<-readLines('AminoAcidNotes.txt')

Amino.Study.A.Notes<-Amino.Study.Notes[grep("Alanine",Amino.Study.Notes)]
Amino.Study.P.Notes<-Amino.Study.Notes[grep("Proline",Amino.Study.Notes)]
Amino.Study.S.Notes<-Amino.Study.Notes[grep("Serine",Amino.Study.Notes)]
Amino.Study.T.Notes<-Amino.Study.Notes[grep("Theronine",Amino.Study.Notes)]

#---------------------------------------------Combinatorics----------------------------------------------
Amino.Study.Peptides.2<-combn(Amino.Study, 2)
Amino.Study.Peptides.3<-combn(Amino.Study, 4)
Amino.Study.Combinations.sequence<-str_c(as.character(Amino.Study.Peptides.2[,2]),collapse="")
Amino.Study.Combinations.sequence.2<-str_c(as.character(Amino.Study.Peptides.3[,2]),collapse="")

Octapeptide<-stri_join(sample(Amino.Study,8,replace=TRUE),collapse="")

Amino.Study.Combinations.sequence.Property.Charge.pH<-lapply(Amino.Study.Combinations.sequence,
                                                             function(x) charge(as.character(x),
                                                                                pH=seq(from = 5,to = 9,by = 0.5), 
                                                                                pKscale="Bjellqvist"))
Amino.Study.Combinations.sequence.Property.Charge.pH.1<-lapply(Amino.Study.Combinations.sequence.2,
                                                               function(x) charge(as.character(x),
                                                                                  pH=seq(from = 5,to = 9,by = 0.5), 
                                                                                  pKscale="Bjellqvist"))
Amino.Study.Combinations.sequence.Property.Charge.pH.3<-lapply(Octapeptide,
                                                               function(x) charge(as.character(x),
                                                                                  pH=seq(from = 5,to = 9,by = 0.5), 
                                                                                  pKscale="Bjellqvist"))

Amino.Study.Combinations.sequence.Property.membrame<-membpos(Amino.Study.Combinations.sequence)
#---------------------------------------------Permutations------------------------------------------------
Amino.Study.Permutations.Matrix<-allPerms(Amino.Study)
#---------------------------------------------Temperature, Pressure and pH--------------------------------
temperature <-seq(0,320,2)
pressure<-seq(1,161,1)
pH<-seq(1,17, 0.1)
aa <- seq2aa("Octapeptide",Octapeptide)
ionize.sequence<-ionize.aa(aa, T=25, pH=pH)
#---------------------------------------------Properties--------------------------------------------------

peptide.study.property<-lapply(Octapeptide, function(x) aacomp(x))

sequence.octapeptide.boman<-lapply(Octapeptide, function(x) boman(x))
sequence.octapeptide.instaIndex<-lapply(Octapeptide, function(x) instaIndex(x))
sequence.octapeptide.membpos<-lapply(Octapeptide, function(x) membpos(x,180))
sequence.octapeptide.mw<-lapply(Octapeptide, function(x) mw(x))

sequence.properties.df<-data.frame()
sequence.properties.df<-cbind("Octapeptide",sequence.octapeptide.boman[[1]],
                              sequence.octapeptide.instaIndex[[1]],
                              sequence.octapeptide.mw[[1]])
colnames(sequence.properties.df)<-c("Type","Boman","instaIndex","MW")
#---------------------------------------------Tables------------------------------------------------------

Table.1<-xtable(amino.study.df)
Table.2<-xtable(sequence.properties.df)

#--------------------------------------------Figures-------------------------------------------------------
Figure.1<-plot(Amino.Study.Combinations.sequence.Property.Charge.pH[[1]],type="l",xlab="pH",
               ylab='Charge',main="Peptide Charge by pH level")
Figure.2<-plot(Amino.Study.Combinations.sequence.Property.Charge.pH.3[[1]],type="l",xlab="pH",
               ylab='Charge',main="Peptide Charge by pH level for Octapeptide")
title("Residue Statistics for the Octapeptide")
Figure.3<-plot(pH,ionize.sequence,type="l",xlab="",ylab='',main="")
Figure.4.NetCharge<-plot(pH, ionize.sequence[, 1], type="l", 
                         xlab="pH",col="blue", 
                         ylab="net charge (Z)")
lines(pH, ionize.aa(aa, T=100, pH=pH)[, 1],lty=2, col="red")
lines(pH, ionize.aa(aa, T=150, pH=pH)[, 1],lty=3, col="green")
title(main="Net Charge Octapeptide by pH and Temperature")
legend("topright",legend=c('T=25','T=100','T=150'),
       cex = 0.7,
       col=c('blue','red','green',
             'yellow','purple','orange','violet'),lty=1:8)
#-------------------------------------------Function Library----------------------------------------------
