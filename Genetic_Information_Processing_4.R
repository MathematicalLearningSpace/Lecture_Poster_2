#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------------------R API -----------------------------------------------
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(rentrez);library(xtable)
library(stringi);library(readr);library(Matrix);library(igraph);library(visNetwork)

#---------------------Data Sets for the Classroom-------------------------------------------------------------------------
#----------------------Formatted Student Notes for the Classroom Lectures-------------------------------------------------
KEGG.IDs.df <- as.data.frame(read_csv("data/KEGG_IDs_For_Query.txt"))
KEGG.Cancer.IDs.df<- as.data.frame(read_csv("data/KEGG_Cancer_IDs_For_Query.txt"))
Genetic.Information.Processing.df <- as.data.frame(read_csv("Data/Genetic_Information_Processing.txt"))
GIP.Folding.Sorting.Degradation.df<-Genetic.Information.Processing.df[9:15,]

#------------------------------------Nucleotide Excision Repair (NER)----------------------------------------------------------

KEGG.Ubiquitin.System.2<-keggGet("path:pxb03420")
KEGG.Ubiquitin.System.2[[1]]$DESCRIPTION

NER.Holo.TFIIH.complex<-c(CAK.subcomplex,TFIIH.core.complex)
NER.CAK.subcomplex<-c("CDK7","XPB","XPD","TTDA","MNATI","CCNH")
NER.TFIIH.core.complex<-c("TFIIH.1","TFIIH.2","TFIIH.3","TFIIH.4")

NER.Cul4.DDB.Complex<-c("RBX1","Cul4","DDB1","DDB2")
NER.XPC.Complex<-c("XPC","HR23B","CETN2")
NER.Cul4.CSA.Complex<-c("RBX1","Cul4","DDB1","CSA","CSB")

NER.Eurkaryotic.Repair.Incision<-c("XPF","ERCC1")
NER.Eurkaryotic.Repair.Excision<-c("Po1.Delta","Po1.epsilon","PCNA","RFC")
NER.Eurkaryotic.Repair.Ligation<-c("Lig1")

NER.GGR.df<-as.data.frame(c(NER.Cul4.DDB.Complex,
           NER.XPC.Complex,
           NER.TFIIH.core.complex,
           NER.Eurkaryotic.Repair.Incision,
           NER.Eurkaryotic.Repair.Excision,
           NER.Eurkaryotic.Repair.Ligation
           ))

NER.TCR.df<-as.data.frame(c(NER.Cul4.CSA.Complex,
           NER.TFIIH.core.complex,
           NER.Eurkaryotic.Repair.Incision,
           NER.Eurkaryotic.Repair.Excision,
           NER.Eurkaryotic.Repair.Ligation
           ))
Nucleotide.Excision.Repair.publications<-entrez_search(db="pubmed", term="Nucleotide Excision Repair", retmax=40)
Holo.TFIIH.complex.publications<-entrez_search(db="pubmed", term="Holo TFIIH complex.", retmax=40)
#----------------------------------------------------------Algorithm----------------------------------------------------------
NER.States<-c("DNA.Lesion",
              "Damage.Recognition",
              "UMP",
              "DNA.Unwinding",
              "Incision",
              "Excision.DNA.Synthesis",
              "Ligation")
#-----------------------------------Mismatch Repair (MR)-----------------------------------------------------------------------
MR.States<-c("Recognition","Excision","Resynthesis","Ligation")

MR.Eurkaryotic.Repair.Recognition<-c("MutL.alpha.PMS2",
                                     "MutL.alpha.MLH1",
                                     "MutS.alpha.MSH6",
                                     "MutS.alpha.MSH2",
                                     "MutL.alpha.MLH1",
                                     "MutL.alpha.PMS2",
                                     "MutS.beta.MSH2",
                                     "MutS.beta.MSH3",
                                     "MLH3.MLH1",
                                     "MLH3.MLH3",
                                     "MutS.beta.MSH2",
                                     "MutS.beta.MSH3")

MR.Eurkaryotic.Repair.Excision<-c("Exol","MutS.homolog","MutL.homolog")
MR.Eurkaryotic.Repair.DNA.resynthesis<-c("RPA","Pol.delta")
MR.Eurkaryotic.Repair.Ligation<-c("LigI")

Nucleotide.Excision.Repair.MutL.alpha.PMS2.publications<-entrez_search(db="pubmed", term="MutL alpha PMS2", retmax=40)
#-----------------------------------------------Algorithm-------------------------------------------

MR.Eurkaryotic.Repair.df<-as.data.frame(c(MR.Eurkaryotic.Repair.Recognition,
                         MR.Eurkaryotic.Repair.Excision,
                         MR.Eurkaryotic.Repair.DNA.resynthesis,
                         MR.Eurkaryotic.Repair.Ligation))


#-------------------------------Homologous.Recombination (HR)------------------------------------------------------------

HR.Eurkarytoic.Strand.Double.Break.Repair.Recognition<-c("ATM","BRCA1",
                                                         "TOPBP1","BRIP1","BARD1","CtIP",MRN.Complex,"Abraxas","RAP80","NBA1",
                                                         "BRE","BRCC36","PALB","DSS1","BRCA2","SYCP3")
Rad51.paralogs<-c("Rad.51.B","Rad.51.C","Rad.51.D","XRCC2","XRCC3")
RAD51.paralogs.publications<-entrez_search(db="pubmed", term="Rad51 paralogs", retmax=40)

MRN.complex<-c("Rad50","Mre11","Nbs1")

HR.Eurkarytoic.Strand.Double.Break.Repair.Formation.Filament<-c("RPA","Rad51","Rad52")
HR.Eurkarytoic.Strand.Double.Break.Repair.Strand.Invasion<-c("RAD54")
HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.SDSA<-c("pol.delta")
HR.Eurkarytoic.Strand.Double.Break.Repair.Strand.displacement<-c("")
HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.SDSA.ligation<-c("BLM")
HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.DSBR<-c("BLM","Mus81","TOP3","Eme1")
HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.BIR<-c("")

Homologous.Recombination.publications<-entrez_search(db="pubmed", term="Homologous.Recombination", retmax=40)

#---------------------------Algorithm---------------------------------------------------------------------------------------

HR.states<-c("Recognition","Filment","Invasion","Synthesis.DSBR","Synthesis.SDSA","Synthesis.BIR")

HR.Eurkarytoic.Strand.Double.Break.Repair.df<-as.data.frame(c(
  HR.Eurkarytoic.Strand.Double.Break.Repair.Recognition,
  HR.Eurkarytoic.Strand.Double.Break.Repair.Formation.Filament,
  HR.Eurkarytoic.Strand.Double.Break.Repair.Strand.Invasion,
  HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.DSBR,
  HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.SDSA,
  HR.Eurkarytoic.Strand.Double.Break.Repair.DNA.Synthesis.BIR))

Eurkarytoic.Strand.Double.Break.Repair.Formation.Filament.publications<-entrez_search(db="pubmed", term="Eurkarytoic.Strand.Double.Break.Repair.Formation.Filament", retmax=40)

#-----------------------------------------------Tables-------------------------------------------------------------------------
Table.1<-xtable(Genetic.Information.Processing.df)

Table.2<-HR.Eurkarytoic.Strand.Double.Break.Repair.df

#-----------------------Figures to be Added by Students-------------------------------------------------------------

Figure.1<-plot()

#-----------------------------------------------References-----------------------------------
Reference.1<-c("","","")
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
