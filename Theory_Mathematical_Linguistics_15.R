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
