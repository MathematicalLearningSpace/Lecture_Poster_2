#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-------------------------------------------R API --------------------------------------------------
library(seqinr);library(ape);library(xtable);library(Biostrings);library(protr)
library(Peptides);library(HDMD);library(genoPlotR);library(SetRank);library(stringr);library(stringi)

detach("package:Biostrings", unload=TRUE)
#---------------Generate Theoretical Data Set Based on Alphabet--------------------------------------------------------
nucleotides <- c("A", "C", "G", "T")
probabilities.1 <-c(0.2, 0.3, 0.3, 0.2) 
probabilities.2<- c(0.2, 0.3, 0.3, 0.2)
probabilities.3<- c(0.2, 0.3, 0.3, 0.2)
probabilities.4<- c(0.2, 0.3, 0.3, 0.2)

seqlength <- 100
DNA.Sequence.1<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.1)
DNA.Sequence.2<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.2)
DNA.Sequence.3<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.3)
DNA.Sequence.4<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.4)
#----------------------Translate Data to Amino Acids Data Type----------------------------------------------------

AA.DNA.Sequence.1<-c2s(translate(as.character(DNA.Sequence.1),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.2<-c2s(translate(as.character(DNA.Sequence.2),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.3<-c2s(translate(as.character(DNA.Sequence.3),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.4<-c2s(translate(as.character(DNA.Sequence.4),frame=0,sens='F',numcode=1))

Gene.Study.Names<-c("AA.DNA.Sequence.1","AA.DNA.Sequence.2","AA.DNA.Sequence.3","AA.DNA.Sequence.4")
#---------Example Residue Properties for Classroom Discussion, Function Modification and Table Design------------------------
Index.Boman<-c(boman(AA.DNA.Sequence.1))
Index.Stability<-c(instaIndex(AA.DNA.Sequence.1))
Index.Aliphatic<-c(aIndex(AA.DNA.Sequence.1))

Gene.Index.Table.df<-data.frame()
Gene.Index.Table.df<-cbind(Index.Boman,Index.Stability,Index.Aliphatic)
colnames(Gene.Index.Table.df)<-c("Boman","Stability","Aliphatic")
#-----------------------------------Tables------------------------------------------------------
Table.1<-xtable(Gene.Index.Table.df)
Table.2<-table(s2c(AA.DNA.Sequence.1)) -> AA.DNA.Sequence.1.tmp
names(AA.DNA.Sequence.1) <- aaa(names(AA.DNA.Sequence.1))

#----------Figures for Classroom Presentation------------------------------------------------------

Figure.1<-dotchart(AA.DNA.Sequence.1.tmp,pch=19,main='Frequency of Amino Acids', xlab='Frequency', ylab='Amino Acid')
Figure.2<-slidingwindowplot(1,s2c(AA.DNA.Sequence.1))

#---------------------Function Library for Modification in Classroom---------------------------------------------
slidingwindowplot <- function(windowsize, inputseq) 
  {                 starts <- seq(1, length(inputseq)-windowsize, by = windowsize) 
                    n <- length(starts) 
                    chunkGCs <- numeric(n) 
                    for (i in 1:n) 
                    { chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
                      chunkGC <- GC(chunk) 
                      chunkGCs[i] <- chunkGC } 
                      # Additional Figures and Tables to be placed here with caption
                      plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}
#-------------Function Template Library for Classroom Presentation and Modification---------------------
f.1<-function(X)
 {
  Z<-""
  a<-1
  W<-runif(length(X),0,1)
  for(i in 1:length(X))
  {  
	Z<-stringr::str_c(Z,X[i])
	W[i]<-a*W[i]
  }
  output<-list()
  output$X<-X
  output$a<-a
  output$Z<-Z
  output$W<-W
  return(output)
 } 
test.f.1<-f.1(letters)
test.f.1
