#----------------------------------Example for Students in the Classroom----------------------------
library(easyPubMed);library(bio3d);library(readr);library(CHNOSZ);library(stringr);library(Peptides);library(Biostrings)
library(seqinr);library(seqLogo);library(msa);library(ape);library(dtw);library(dtwclust);library(odseq);library(rphast)
library(plyr)
#-----------------------------------------Data--------------------------------------------------------
data(thermo)
data(aaindex)
#-----------------------------------------Search Engine Design----------------------------------------
Proteasome_xml <- fetch_pubmed_data(get_pubmed_ids("Proteasome"))
Chaperone_xml <- fetch_pubmed_data(get_pubmed_ids("Chaperone"))
Ribosome_xml <- fetch_pubmed_data(get_pubmed_ids("Ribosome"))
Heat_Shock_Proteins<- fetch_pubmed_data(get_pubmed_ids("Heat Shock Proteins"))
Ligand_xml<- fetch_pubmed_data(get_pubmed_ids("Ligand"))
Entropy_xml<- fetch_pubmed_data(get_pubmed_ids("Entropy"))
Lysosome_xml<-fetch_pubmed_data(get_pubmed_ids("Lysosome"))
Endosome_xml<-fetch_pubmed_data(get_pubmed_ids("Endosome"))
Enzymes.Proteins.lysosome_xml<-fetch_pubmed_data(get_pubmed_ids("Lysosome Enzymes Proteins"))
enzyme.inhibitors_xml<-fetch_pubmed_data(get_pubmed_ids("enzyme inhibitors"))
enzyme.active.site_xml<-fetch_pubmed_data(get_pubmed_ids("three-dimensional structure enzyme active site"))
#----------------------------------------------------------------------------------------------------------
ADP.glucose.pyrophosphorylase<-fetch_pubmed_data(get_pubmed_ids("ADP-glucose pyrophosphorylase"))
peptidyl.tyrosine.autophosphorylation_xml<-fetch_pubmed_data(get_pubmed_ids("peptidyl-tyrosine autophosphorylation"))
positive.regulation.of.calcium.mediated.signaling<-fetch_pubmed_data(get_pubmed_ids("positive regulation of calcium mediated signaling"))
mRNA.alternative.polyadenylation<-fetch_pubmed_data(get_pubmed_ids("mRNA alternative polyadenylation"))
Delta.3.UTR<-fetch_pubmed_data(get_pubmed_ids("profiling the changes of 3′-UTR length"))
ER.Stress<-fetch_pubmed_data(get_pubmed_ids("Endoplasmic reticulum stress protein glycosylation"))
InMAP<-fetch_pubmed_data(get_pubmed_ids("Integrative Model for Alternative Polyadenylation"))
mRNA.post.transcriptional.regulation<-fetch_pubmed_data(get_pubmed_ids("mRNA post-transcriptional regulation")) 
Notes.content.A(mRNA.post.transcriptional.regulation,"mRNA post transcriptional regulation","cancer",FALSE,save.notes=TRUE)
#--------------------------------------------------------------------------------------------------------
#--------------------------------------Signal Transduction Networks---------------------------------------
A.xml<-fetch_pubmed_data(get_pubmed_ids("Ras signaling pathway")) 
B.xml<-fetch_pubmed_data(get_pubmed_ids("PI3K-Akt signaling pathway"))  
C.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway")) 
D.xml<-fetch_pubmed_data(get_pubmed_ids("Wnt signaling pathway")) 
E.xml<-fetch_pubmed_data(get_pubmed_ids("p53 signaling pathway"))  
F.xml<-fetch_pubmed_data(get_pubmed_ids("Thyroid hormone signaling pathway")) 
G.xml<-fetch_pubmed_data(get_pubmed_ids("FoxO signaling pathway"))
H.xml<-fetch_pubmed_data(get_pubmed_ids("mTOR signaling pathway")) 
I.xml<-fetch_pubmed_data(get_pubmed_ids("ErbB signaling pathway"))  
J.xml<-fetch_pubmed_data(get_pubmed_ids("MAPK signaling pathway"))  
K.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway"))  
L.xml<-fetch_pubmed_data(get_pubmed_ids("Apoptosis"))  
M.xml<-fetch_pubmed_data(get_pubmed_ids("Cell cycle")) 
N.xml<-fetch_pubmed_data(get_pubmed_ids("Rap1 signaling pathway"))  
O.xml<-fetch_pubmed_data(get_pubmed_ids("Chemokine signaling pathway"))
p.xml<-fetch_pubmed_data(get_pubmed_ids("Phosphatidylinositol signaling system"))
Q.xml<-fetch_pubmed_data(get_pubmed_ids("T cell receptor signaling pathway"))
R.xml<-fetch_pubmed_data(get_pubmed_ids("Estrogen signaling pathway"))
S.xml<-fetch_pubmed_data(get_pubmed_ids("VEGF signaling pathway"))
genes.signaling.pathway<-c("NF1","NRAS", "PIK3CA", "NRAS", "PTEN", "TP53", "CTNNB1", "SMAD4", 
                           "CDKN2A", "PTEN", "TP53", "NCOR1", "BRAF", "SMAD4")
for(i in 1: length(genes.signaling.pathway))
{
genes.xml<-fetch_pubmed_data(get_pubmed_ids(genes.signaling.pathway[i]))   
Notes.content.A(genes.xml,stringr::str_c("Genes signaling_",genes.signaling.pathway[i]),"cancer",FALSE,save.notes=TRUE)
}
Notes.content.A(A.xml,"Ras signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(B.xml,"PI3K-Akt signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(C.xml,"Hippo signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(D.xml,"Wnt signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(E.xml,"p53 signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(F.xml,"Thyroid hormone signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(G.xml,"FoxO signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(H.xml,"mTOR signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(I.xml,"ErbB signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(J.xml,"MAPK signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(K.xml,"Hippo signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(L.xml,"Apoptosis","cancer",FALSE,save.notes=TRUE)
Notes.content.A(M.xml,"Cell cycle","cancer",FALSE,save.notes=TRUE)
Notes.content.A(N.xml,"Rap1 signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(O.xml,"Chemokine signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(P.xml,"Phosphatidylinositol signaling system","cancer",FALSE,save.notes=TRUE)
Notes.content.A(Q.xml,"T cell receptor signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(R.xml,"Estrogen signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(S.xml,"VEGF signaling pathway","cancer",FALSE,save.notes=TRUE)
Notes.content.A(genes.xml,"Genes signaling","cancer",FALSE,save.notes=TRUE)
#-------------------------------------------------------------------------------------------------------------------
Cell.based.Mathematical.Model_xml<-fetch_pubmed_data(get_pubmed_ids("Cell based Mathematical Model Protein Complex"))
molecular.docking_xml<-fetch_pubmed_data(get_pubmed_ids("molecular docking"))
molecular.machine_xml<-fetch_pubmed_data(get_pubmed_ids("molecular machine"))
molecular.mechanics_xml<-fetch_pubmed_data(get_pubmed_ids("molecular mechanics"))
Reversible.phosphorylation.conformational.change_xml<- fetch_pubmed_data(get_pubmed_ids("Reversible phosphorylation conformational change"))
Crizotinib_xml<-fetch_pubmed_data(get_pubmed_ids("Crizotinib"))
Pyrimidine.dimers_xml<-fetch_pubmed_data(get_pubmed_ids("Pyrimidine dimers"))
intracellular.pH.temperature_xml<-fetch_pubmed_data(get_pubmed_ids("intracellular pH temperature"))
c.Met.inhibitors_xml<-fetch_pubmed_data(get_pubmed_ids("C-Met Inhibitors"))
#-----------------------------G Quadruplexes-------------------------------------------
Intramolecular.G.Quadruplexes.Basket_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Basket"))
Intramolecular.G.Quadruplexes.Chair_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Chair"))
Intramolecular.G.Quadruplexes.Propeller_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Propeller"))
Intermolecular.G.Quadruplexes.Hairpin.dimmer_xml<-fetch_pubmed_data(get_pubmed_ids("Intermolecular G-Quadruplexes Hairpin dimmer"))
Intermolecular.G.Quadruplexes.Tetrameric_xml<-fetch_pubmed_data(get_pubmed_ids("Intermolecular G-Quadruplexes Tetrameric")
#-----------------------Mathematical Models--------------------------------------------
Differential.Equations_xml<-fetch_pubmed_data(get_pubmed_ids("Differential Equations"))
MathModels_xml<-fetch_pubmed_data(get_pubmed_ids("Mathematical Models"))
                                                                
#-----------------------------------------------------------------------Reading Lists----------------------------------------------
reading.list(Lysosome_xml,"cancer",TRUE,"Lysosome")
reading.list(Enzymes.Proteins.lysosome_xml,"cancer",TRUE,"Lysosome_Enzymes_Proteins")
reading.list(MathModels_xml,"cancer",TRUE,"MathModels")
reading.list(Cell.based.Mathematical.Model_xml,"cancer",TRUE,"Cell_based_Mathematical_Model")
reading.list(enzyme.inhibitors_xml,"cancer",TRUE,"enzyme_inhibitors")
reading.list(enzyme.active.site_xml,"cancer",TRUE,"three_dimensional_structure_enzyme_active_site")
reading.list(molecular.docking_xml,"cancer",TRUE,"molecular docking")
reading.list(molecular.machine_xml,"cancer",TRUE,"molecular machine")
reading.list(molecular.mechanics_xml,"cancer",TRUE,"molecular mechanics")
reading.list(Reversible.phosphorylation.conformational.change_xml,"cancer",TRUE,"Reversible_phosphorylation_conformational_change")
reading.list(intracellular.pH.temperature_xml,"cancer",TRUE,"intracellular_pH_temperature")
reading.list(Pyrimidine.dimers_xml,"cancer",TRUE,"Pyrimidine_dimers")
reading.list(Crizotinib_xml,"cancer",TRUE,"Crizotinib")
#------------------------------Notes Content--------------------------------------------------
Notes.content.A(Enzymes.Proteins.lysosome_xml,"Enzymes_Proteins_Lysosome","cancer",save.notes=TRUE)
Notes.content.A(MathModels_xml,"MathModels","cancer",save.notes=TRUE)
Notes.content.A(Pyrimidine.dimers_xml,"Dimers","cancer",save.notes=TRUE)
Notes.content.A(molecular.machine_xml,"Molecular_Machine","cancer",save.notes=TRUE)
Notes.content.A(molecular.docking_xml,"Molecular_Docking","cancer",save.notes=TRUE)
Notes.content.A(molecular.mechanics_xml,"Molecular_Mechanics","cancer",save.notes=TRUE)
Notes.content.A(Endosome_xml,"Endosome","cancer",save.notes=TRUE)
Notes.content.A(Crizotinib_xml,"Crizotinib","cancer",save.notes=TRUE)
Notes.content.A(c.Met.inhibitors_xml,"C_Met Inhibitors","cancer",save.notes=TRUE)
#---------------------------------------------------------------------------------------------
Notes.content.A(Intramolecular.G.Quadruplexes.Basket_xml,"Intramolecular_G_Quadruplexes_Basket","cancer",FALSE,save.notes=TRUE)
Notes.content.A(Intramolecular.G.Quadruplexes.Chair_xml,"Intramolecular_G_Quadruplexes_Chair","cancer",FALSE,save.notes=TRUE)
Notes.content.A(Intramolecular.G.Quadruplexes.Propeller_xml,"Intramolecular_G_Quadruplexes_Propeller","cancer",FALSE,save.notes=TRUE)
Notes.content.A(Intermolecular.G.Quadruplexes.Hairpin.dimmer_xml,"Intermolecular_G_Quadruplexes_Hairpin_dimmer","cancer",FALSE,save.notes=TRUE)
Notes.content.A(Intermolecular.G.Quadruplexes.Tetrameric_xml,"Intermolecular_G_Quadruplexes_Tetrameric","cancer",FALSE,save.notes=TRUE)
                                                                
Notes.content.A(Differential.Equations_xml,"Differential_Equations","cancer",FALSE,save.notes=TRUE)

                                                         

#------------------------------------Table Designs----------------------------------------------------
#------------------Ribosome, Chaperone and Proteasome Ligands from RCSB PDB-----------
Table.1 <- read_csv("Ligands.csv")
#------------------------------------Figure Designs--------------------------------------------------


#------------------------------------Function Library Designs----------------------------------------
reading.list<-function(X,search.terms,save.notes=FALSE, notes.name)
{
 X_titles <- unlist(xpathApply(X, "//ArticleTitle", saveXML))
  X_titles_1 <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", X_titles)
  #X_titles_1[nchar(X_titles_1)>75] <- paste(substr(X_titles_1[nchar(X_titles_1)>75], 1, 70), "...", sep = "")
  if(rlfilter){X_title_1_filtered<-X_titles_1[grep(search.terms,X_titles_1)]}
  if(save.notes)
  {
    write(X_title_1_filtered, file=stringr::str_c("Notes_Reading_List_",notes.name,"_",search.terms,".txt"),append=FALSE)
  }
  output<-list()
  output$Titles<-X_titles_1
  output$Titles.Filtered<-X_title_1_filtered
  return(output)
}
reading.list(Proteasome_xml,"ligand",TRUE,"Proteasome")
reading.list(Chaperone_xml,"ligand",TRUE,"Chaperone") 
reading.list(Ribosome_xml,"ligand",TRUE,"Ribosome") 
reading.list(Heat_Shock_Proteins,"ligand",TRUE,"Heat_Shock_Proteins")
reading.list(Ligand_xml,"Receptor Tryosine Kinases",TRUE,"Ligand")
reading.list(Entropy_xml,"ligand",TRUE,"Entropy")

#----------------------------------------------------------Notes Content-------------------------------------------
Notes.content.A<-function(X,note.name,search.terms,rlfilter=TRUE,save.notes=FALSE)
{
  X.list<-list()
  X_abstract <- unlist(xpathApply(X, "//AbstractText", saveXML))
  X_abstract_1 <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", X_abstract)  
  X_titles <- unlist(xpathApply(X, "//ArticleTitle", saveXML))
  X_titles_1 <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", X_titles)
  if(rlfilter){
    X_abstract_1_filtered<-  X_abstract_1[grep(search.terms,X_abstract_1)]
    X_titles_1_filtered<-X_titles_1[grep(search.terms,X_titles_1)]
  }
  else{X_abstract_1_filtered<-X_abstract_1
  X_titles_1_filtered<-X_titles_1}
  notes.1<-X_abstract_1_filtered
  readingList<-X_titles_1_filtered
  notes.summary.format<-""
  readingList.format<-""
  for(i in 1:length(X_abstract_1_filtered))
  {
    X.list[[i]]<-strsplit(notes.1[i],". ",fixed=TRUE)
    notes.summary<-unlist(X.list[[i]])
    for(j in 1:length(notes.summary))
    {
      notes.summary.format<-stringr::str_c(notes.summary.format, "\\item ", notes.summary[j], " \\cite{key",i,"}   " )
    }
    notes.summary.format<-stringr::str_c(notes.summary.format,"%---------------------------------------------------%")
    
  }
  for(i in 1:length(readingList))
  {
    readingList.format<-stringr::str_c(readingList.format," \\bibitem[",i,"]{key",i,"}", readingList[i])
  }
  notes.and.reading.list<-stringr::str_c(notes.summary.format,"%-------------References--------------% ",readingList.format)
  if(save.notes)
  {
    write(notes.and.reading.list, file=stringr::str_c("Classroom_Notes_For_Journal_Articles_",note.name,"_",search.terms,".txt"),append=FALSE)
  }
  
  output<-list()
  output$Notes<-notes.1
  output$Notes.Summary<-notes.summary.format
  output$Reading.List<-readingList.format
  return(output)
}

test.Notes.content.A<-Notes.content.A(Lysosome_xml,"Lysosome","cancer",save.notes=TRUE)
test.Notes.content.A

#------------------------------------References------------------------------------------------------

#H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
#(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
