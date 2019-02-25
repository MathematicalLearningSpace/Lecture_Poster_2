#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------------------------------R API ------------------------------------------------
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(rentrez)
library(xtable);library(stringi);library(readr);library(readxl);library(Matrix);library(igraph);library(visNetwork)
#-------------------------------------------------Data-------------------------------------------------------------------------
#--------------Design Data Tables to be used in the Classroom------------------------------------------------------------------
KEGG.IDs.df <- as.data.frame(read_csv("data/KEGG_IDs_For_Query.txt"))
KEGG.Cancer.IDs.df<- as.data.frame(read_csv("data/KEGG_Cancer_IDs_For_Query.txt"))

View(KEGG.IDs.df)

Genetic.Information.Processing.df <- as.data.frame(read_csv("Data/Genetic_Information_Processing.txt"))
View(Genetic.Information.Processing.df)

KEGG.organisms<-keggList("organism")
KEGG.Brite<-keggList("brite")
KEGG.organisms.df<-as.data.frame(KEGG.organisms)
KEGG.organisms.AT<-KEGG.organisms.df[177,]

EUtils.Database <- entrez_info()
PMC.database<-entrez_info(db = "pmc")

uniprot.species.human.df <- as.data.frame(read_excel("uniprot-all.xlsx"))
View(uniprot.species.human.df)

uniprot.species.human.CDK7.df<-uniprot.species.human.df[grep("CDK7",uniprot.species.human.df$`Gene names`),]
CDK7.pdb<-getUniProt(uniprot.species.human.CDK7.df$Entry[1])

Protein_Report_RCSB_PDB_Results.df <- as.data.frame(read_csv("data/Protein_Report_RCSB_PDB_Results.csv"))
View(Protein_Report_RCSB_PDB_Results.df)
#------------------------------------------------Notes on Ubiquitin System and Proteasome Assembly----------------------------------'
Ubiquitin.System.Notes<-as.data.frame(read_csv("Data/Ubiquitin_System_Notes.txt"))
Proteasome.Assembly.Notes<-as.data.frame(read_csv("Data/Proteasome_Assembly_Notes.txt"))

#----------------------------------------------Ubiquitin System: Ubiquitin mediated proteolysis---------------------------------------------------------------
#
KEGG.Ubiquitin.System<-keggGet("br:ko04121")
KEGG.Ubiquitin.System.1<-keggGet("path:pxb04120")
KEGG.Ubiquitin.System.1[[1]]$DESCRIPTION

KEGG.Ubiquitin.System.1.Gene<-KEGG.Ubiquitin.System.1[[1]]$GENE
#
#---------------------------------------------B: Standard 20s Proteasome AT----------------------------------

KEGG.Proteasome.20s<-keggGet("path:pxb03050")

KEGG.Proteasome.20s.Description<-KEGG.Proteasome.20s[[1]]$DESCRIPTION
KEGG.Proteasome.20s.Gene<-KEGG.Proteasome.20s[[1]]$GENE
KEGG.Proteasome.20s.Articles<-KEGG.Proteasome.20s[[1]]$REFERENCE

number.even.1<-seq(0,length(KEGG.Proteasome.20s[[1]]$GENE),2)
KEGG.Proteasome.20s.Gene<-KEGG.Proteasome.20s[[1]]$GENE[number.even.1]

KEGG.Proteasome.20s.Articles.References<-NULL
for(i in 1:length(KEGG.Proteasome.20s.Articles))
{
  KEGG.Proteasome.20s.Articles.References[i]<-stri_join(KEGG.Proteasome.20s.Articles[[i]]$REFERENCE,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$AUTHORS,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$TITLE,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$JOURNAL)
}
KEGG.Proteasome.20s.Articles.References.df<-as.data.frame(KEGG.Proteasome.20s.Articles.References)
colnames(KEGG.Proteasome.20s.Articles.References.df)<-c("Reference")

Proteasome.20s.df<-data.frame()
Protesome.core.particle.alpha<-c("alpha.1","alpha.2","alpha.3","alpha.4","alpha.5","alpha.6","alpha.7")
Protesome.core.particle.beta<-c("beta.1","beta.2","beta.3","beta.4","beta.5","beta.6","beta.7")
#Reverse the sequence
Protesome.core.particle.alpha.reverse<-c("alpha.7","alpha.6","alpha.5","alpha.4","alpha.3","alpha.6","alpha.7")
Protesome.core.particle.beta.reverse<-c("beta.7","beta.6","beta.5","beta.4","beta.3","beta.2","beta.1")
Proteasome.20s.Middle<-c(Protesome.core.particle.alpha,
                         Protesome.core.particle.beta,
                         Protesome.core.particle.beta.reverse,
                         Protesome.core.particle.alpha.reverse)
Proteasome.20s.Assembled.df<-as.data.frame(c(Proteasome.20s.Middle))

#----------------------------------------------A: PA 700 26s Proteasome Assembly------------------------------------------------------------

KEGG.Proteasome.26s.RPN.13.publications<-entrez_search(db="pubmed", term="26S proteasome regulatory subunit RPN13-like", retmax=40)
Proteasome.26s.df<-data.frame()
Proteasome.26s.lid<-c("RPN.3","RPN.5","RPN.6","RPN.7","RPN.8","RPN.9","RPN.11","RPN.12","RPN.15")
Proteasome.26s.base<-c("RPN.1","RPN.2","RPN.13","RPT.1","RPT.2","RPT.6","RPT.4","RPT.5","RPT.3")
Proteasome.26s.Top<-c(Proteasome.26s.lid,Proteasome.26s.base)
Proteasome.26s.Bottom<-c(Proteasome.26s.base,Proteasome.26s.lid)
Proteasome.26s.Assembled.df<-as.data.frame(c(Proteasome.26s.Top,
                                             Proteasome.20s.Middle,
                                             Proteasome.26s.Bottom))
#--------------------------------------------C: Immunoproteasome-------------------------------------------

Immunoproteasome.20s.df<-data.frame()
Immunoproteasome.20s.core.Regulatory.Particle<-c("PA28.alpha","PA28.beta") # hetro hexamer or heptamer
Immunoproteasome.20s.core.Regulatory.Particle.1<-c("PA28.gamma") #, homo hexamer
Immunoproteasome.20s.core.particle.beta<-c("beta.5i","beta.2i","beta.1i")
Immunoproteasome.20s.core.Regulatory.Particle.reverse<-c("PA28.alpha","PA28.beta")
Immunoproteasome.20s.core.particle.beta.reverse<-c("beta.1i","beta.2i","beta.5i")

Immunoproteasome.20s.Assembled.df<-as.data.frame(c(Immunoproteasome.20s.core.Regulatory.Particle,
                                  Immunoproteasome.20s.core.particle.beta,
                                  Proteasome.20s.Middle,
                                  Immunoproteasome.20s.core.particle.beta.reverse,
                                  Immunoproteasome.20s.core.Regulatory.Particle.reverse))

Thymoproteasome.df<-data.frame()
Thymoproteasome.core.particle.beta<-c("beta.5t")
Thymoproteasome.Assembled.df<-as.data.frame(c(Proteasome.26s.Top,Proteasome.20s.Middle,Proteasome.26s.Bottom))

#-----------------------------------------------Tables-------------------------------------------------------------------------

Table.1<-xtable(Genetic.Information.Processing.df)
Table.2<-xtable(Proteasome.26s.Assembled.df)
Table.3<-xtable(Proteasome.20s.Assembled.df)
Table.4<-xtable(Immunoproteasome.20s.Assembled.df)
Table.5<-xtable(Thymoproteasome.Assembled.df)

#-------------Figures to be designed by students in the classroom--------------------------------------------------------------

Figure.1<-plot()
#-----------------------------------------------References-----------------------------------
References<-xtable(KEGG.Proteasome.20s.Articles.References.df)
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



