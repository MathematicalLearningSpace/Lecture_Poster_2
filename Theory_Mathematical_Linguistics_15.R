#----------------------------------Example for Students in the Classroom----------------------------
library(easyPubMed)
library(bio3d)
library(readr)

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

#------------------------------------References------------------------------------------------------

#H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
#(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
