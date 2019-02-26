#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-------------------------------------------R API ----------------------------------------------------------
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(rentrez);library(xtable)
library(stringi);library(readr);library(Matrix);library(igraph);library(visNetwork);library(Hmm);library(markovchain)

#-------------------Data Sets for the Classroom-------------------------------------------------------------------------
#-------------------Formatted Student Notes for the Classroom Lectures--------------------------------------------------
KEGG.IDs.df <- as.data.frame(read_csv("data/KEGG_IDs_For_Query.txt"))
KEGG.Cancer.IDs.df<- as.data.frame(read_csv("data/KEGG_Cancer_IDs_For_Query.txt"))
Genetic.Information.Processing.df <- as.data.frame(read_csv("Data/Genetic_Information_Processing.txt"))
GIP.Folding.Sorting.Degradation.df<-Genetic.Information.Processing.df[9:15,]

Protein.synthesis.attenuating<-c("")

#-----------------------------------------------Protein Processing in Endoplasmic Reticulum------------------------------------------

Ribosome.Protein.Export<-c("p180","OSTs","G1cI","Climp63","Sec61","Sec62/63")

ER.unfolded.lectin<-c("G1cI","G1cII","CNX","ERP57","CRT")
ER.folded.intermediate<-c("G1M9")
ER.folded<-c("M9","ERManI","VIP36","ERGIC53","ERGL")
ER.folded.Golgi.Cop.2<-c("SAR1","Sec13/31","Sec23/24","Sec12")
ER.misfolded.apoptosis<-c("BAK/BAX","Bcl2","Calpain","CASP12")
ER.misfolded.terminal.target.protein<-c("Erol1","PDIs","OS9","XTP3B","Bap31","TRAP","Sec61")

ER.Stress.UPR<-c("PERK","eIF2.alpha","NRF2","GADD34","ATF4")
ER.Stress.Misfolded<-c("IRE1","TRAF2","ASK1","MKK7","JNK")
ER.Stress.Misfolded.sum<-c("ATF6","COPII","WFS1","ATF6.p50")
ER.Stress.Nucleus<-c("AARE","ERSE","URPE")
ER.Stress.Nucleus.metabolism<-c("CHOP","Bcl2")
                       
ERAD.Ubiquitin.ligase.complex<-c("Ubx","TRAM","Derlin","SVIP","p97","NP14","Ufd1")
ERAD<-c("Hsp70","Hsp40","Hsp90","NEF","sHSF","Otu1","DOA1","DSK2","RAD23","Pngl","Ufd2","DUB")

ER.Ubiquitin.ligase.complex.Cytoplasm<-c("UbcH5","Hsp40","Hsp70","CHIP","Parkin")
ER.Ubiquitin.ligase.complex.membrame.ERAD.C<-c("Doa10","Ubc6/7","Cue1","gp78","RMA1","UBE2G2","Derlin")
ER.Ubiquitin.ligase.complex.ERAD.L.M<-c("VIMP","Sel1L","Derlin1","HERP","HRD1","Ubc6/7")

#-----------------------------------------------Algorithm-------------------------------------------
ER.States<-c("Entry",
             "Unfolded.1",
             "Unfolded.2",
             "Folded",
             "Misfolded",
             "Misfolded.Terminal",
             "Exit")

ERAD.Proteasome<-c()

ER.Protein.Path.1.df<-as.data.frame(c(Ribosome.Protein.Export,
                     ER.unfolded.lectin,
                     ER.folded.intermediate,
                     ER.folded,
                     ER.folded.Golgi.Cop.2,
                     ERAD.Proteasome))

ER.Protein.Path.2.df<-as.data.frame(c(Ribosome.Protein.Export,
                     ER.folded.intermediate,
                     ER.Stress.Misfolded,
                     ER.Stress.Misfolded.sum,
                     ER.Stress.Nucleus,
                     ER.misfolded.terminal.target.protein,
                     ERAD.Proteasome))
#-----------------------------------------------Tables-------------------------------------------------------------------------

Table.1<-xtable(Genetic.Information.Processing.df)
Table.2<-xtable(ER.Protein.Path.1.df)
Table.3<-xtable(ER.Protein.Path.2.df)

#-----------------Figures to be Designed by Students in the Classroom----------------------------------------------------------

Figure.1<-plot()

#-----------------References---------------------------------------------------------------------

#----------------------------------------------Function Library----------------------------------------------------------------
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
