#----------------------------------Example for Students in the Classroom----------------------------
library(easyPubMed)
library(bio3d)
library(readr)
library(CHNOSZ)
library(stringr)
library(Peptides)
library(Biostrings)
library(seqinr)
library(seqLogo)
library(msa)
library(ape)
library(dtw)
library(dtwclust)
library(odseq)
library(rphast)
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

#-----------------------------------List of Article Titles---------------------------------------------
Proteasome_titles <- unlist(xpathApply(Proteasome_xml, "//ArticleTitle", saveXML))
Proteasome_titles_1 <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", Proteasome_titles)
Proteasome_titles_1[nchar(Proteasome_titles_1)>75] <- paste(substr(Proteasome_titles_1[nchar(Proteasome_titles_1)>75], 1, 70), "...", sep = "")
print(Proteasome_titles_1)

#------------------------------------Table Designs----------------------------------------------------
#------------------Ribosome, Chaperone and Proteasome Ligands from RCSB PDB-----------
Table.1 <- read_csv("Ligands.csv")
#------------------------------------Figure Designs--------------------------------------------------


#------------------------------------Function Library Designs----------------------------------------
reading.list<-function(X,search.terms,save.notes=FALSE, notes.name)
{
  X_titles <- unlist(xpathApply(X, "//ArticleTitle", saveXML))
  X_titles_1 <- gsub("(^.{5,10}Title>)|(<\\/.*$)", "", X_titles)
  X_title_1_filtered<-X_titles_1[grep(search.terms,X_titles_1)]
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

#------------------------------------References------------------------------------------------------

#H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
#(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
