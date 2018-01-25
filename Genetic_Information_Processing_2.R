library(KEGG.db)
library(KEGGgraph)
library(KEGGprofile)
library(KEGGREST)
library(rentrez)
library(xtable)
library(stringi)
library(readr)
library(Matrix)
library(igraph)
library(visNetwork)

#-------------------------------------------------Data-------------------------------------------------------------------------
#
KEGG.IDs.df <- as.data.frame(read_csv("data/KEGG_IDs_For_Query.txt"))
KEGG.Cancer.IDs.df<- as.data.frame(read_csv("data/KEGG_Cancer_IDs_For_Query.txt"))
Genetic.Information.Processing.df <- as.data.frame(read_csv("Data/Genetic_Information_Processing.txt"))
GIP.Folding.Sorting.Degradation.df<-Genetic.Information.Processing.df[9:15,]

#-----------------------------------------------Ubiquitin mediated proteolysis------------------------------------------


Ubiquitin.enzyme.activation.E1<-c("UB.E1","UB.L.E1.A","UB.L.E1.B","UB.E1.C")
Ubiquitin.enzyme.conjugation.E2<-c("UB.E2.A","UB.E2.B","UB.E2.C","UB.E2.D",
                                   "UB.E2.E","UB.E2.F","UB.E2.G.1","UB.E2.G.2",
                                   "UB.E2.H","UB.E2.I","UB.E2.J.1","UB.E2.J.2",
                                   "UB.E2.L.3","UB.E2.L.6","UB.E2.M","UB.E2.N",
                                   "UB.E2.O","UB.E2.Q","UB.E2.R","UB.E2.S",
                                   "UB.E2.U","UB.E2.W","UB.E2.Z","HIP2","APOLLON")

Ubiquitin.ligase.E3.HECT<-c("UB.3.B","UB.E.3.C")
Ubiquitin.ligase.E3.U-Box<-c("UB.E4.A","UB.E4.B","CHIP","CYC4","PRP19","UIP5")
Ubiquitin.ligase.E3.RING.Finger.Single<-c("Mdm2","BRCA1","COP1")
Ubiquitin.ligase.E3.RING.Finger.Subunit.Multi<-c("")
Ubiquitin.APC.subunits.other<-c("Apc.1","Apc.3","Apc.4","Apc.5","Apc.6","Apc.7","Apc.8","Apc.9","Apc.10","Apc.12","Apc.13")
#
#-----------------------------------------------Algorithm---------------------------------------------------------------------
#
#Step 1. Ub connect E1. Result:ATP->AMP
#Step 2. Ub enters domain:E2,E3. Result: Ub connect E2
#Step 3. Ub enters domain:target recognition. Result: UB connect Target, Ub replication
#Step 4. Ub Polyubiquitination Result: Ub replication
#Step 5. Ub enters Proteasome Result: Ub separated Ub.s
#Step 6. Ub.s recycled

Polyubiquitination.publications<-entrez_search(db="pubmed", term="Polyubiquitination", retmax=40)
Ubiquitin.mediated.proteolysis.publications<-entrez_search(db="pubmed", term="Ubiquitin mediated proteolysis", retmax=40)

#-----------------------------------------------Tables-------------------------------------------------------------------------


Table.1<-xtable(Genetic.Information.Processing.df)


#-----------------------------------------------Figures------------------------------------------------------------------------

Figure.1<-plot()

#-----------------------------------------------References-----------------------------------
Reference.1<-c("","","")

#----------------------------------------------Function Library----------------------------------------------------------------

