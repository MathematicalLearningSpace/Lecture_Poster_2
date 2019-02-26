#----------------------------------------------R API ------------------------------------------------------
library(bio3d);library(xtable);library(Peptides);library(stringi);library(stringr);library(xtable);library(readr);library(readxl)
#---------------------------------------------Data From String DB and GeneCards---------------------------------------------
#---------------Formatted Student Notes for Lectures------------------------------
WD64.Enrichment.Process <- as.data.frame(read_delim("enrichment.Process.tsv", "\t", escape_double = FALSE, trim_ws = TRUE))

WD64.Enrichment.Process$`#pathway ID`
WD64.Enrichment.Process$`pathway description`

D6CR1_Human.Notes<-as.data.frame(read_delim('D6RCR1_Human.txt', "\t", escape_double = FALSE, trim_ws = TRUE))

WDR64.String.Interactions.Notes<-as.data.frame(read_delim('string_interactions.tsv', "\t", escape_double = FALSE, trim_ws = TRUE))
WDR64.Protein.Annotations.Notes<-as.data.frame(read_delim('string_protein_annotations.tsv', "\t", escape_double = FALSE, trim_ws = TRUE))

S.phase.kinase_associated.protein.1<-WDR64.Protein.Annotations.Notes$annotation[1]
immune.response_regulating.cell.surface.receptor.signaling.pathway<-WD64.Enrichment.Process[4,]
immune.response_regulating.cell.surface.receptor.signaling.pathway.network<-WD64.Enrichment.Process$`matching proteins in your network (labels)`[4]

D6RCR1.Human.pdb<-read.pdb("D6RCR1_HUMAN_1.pdb")

Proteins.Study<-c(immune.response_regulating.cell.surface.receptor.signaling.pathway.network)

#------------------------------Tables-------------------------------------------

Table.1<-xtable(WDR64.Protein.Annotations.Notes)
Table.2<-xtable(WD64.Enrichment.Process)

#----------------------Figures to be added by students------------------------------------------


#---------------------Function Library to be added by students----------------------------------
