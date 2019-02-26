#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------------------R API ------------------------------------------------------------
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(rentrez);library(xtable)
library(stringi);library(readr);library(Matrix);library(igraph);library(visNetwork);library(protr);library(Biostrings)
library(Decipher);library(msa);library(seqLogo);library(seqinr);library(odseq);library(bio3d);library(ape)
library(ips);library(cba)

#------------------------------------------------------Data-----------------------------------------------------------------
Protein_Report_RCSB_PDB_Results.df <- as.data.frame(read_csv("data/Protein_Report_RCSB_PDB_Results.csv"))
View(Protein_Report_RCSB_PDB_Results.df)

#-----------------------------------------------------Classification Categories---------------------------------------------

System.Classification.Categories<-unique(Protein_Report_RCSB_PDB_Results.df$Classification)
System.Classification.Categories.N<-length(System.Classification.Categories)

#-----------------------------------------------------DNA Repair and Immune System-------------------------------------------

System.Classification.Categories.DNA.Repair<-System.Classification.Categories[grep("DNA REPAIR",System.Classification.Categories)]
System.Classification.Categories.System.Immune<-System.Classification.Categories[grep("IMMUNE SYSTEM",System.Classification.Categories)]

Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df <-Protein_Report_RCSB_PDB_Results.df[grep("IMMUNE SYSTEM/HYDROLASE",Protein_Report_RCSB_PDB_Results.df$Classification),]
Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df <-Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df[grep("A",Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df$`Chain ID`),]

Protein_Report_RCSB_PDB_Results.System.Immune.DNA.Repair.Hydrolase.df <-Protein_Report_RCSB_PDB_Results.df[grep(System.Classification.Categories.DNA.Repair[2],
                                                                                                     Protein_Report_RCSB_PDB_Results.df$Classification),]
Protein_Report_RCSB_PDB_Results.System.Immune.DNA.Repair.Hydrolase.df <-Protein_Report_RCSB_PDB_Results.System.Immune.DNA.Repair.Hydrolase.df[grep("A",Protein_Report_RCSB_PDB_Results.System.Immune.DNA.Repair.Hydrolase.df$`Chain ID`),]

#----------------------------------------------------Immune System Hydrolase--------------------------------------------------

System.Immune.Hydrolase.Sequences.df<-data.frame()
System.Immune.Hydrolase.Sequences.df<-cbind(Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df$`PDB ID`,
              Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df$`Structure Title`,
              Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df$`Residue Count`,
              Protein_Report_RCSB_PDB_Results.System.Immune.Hydrolase.df$Sequence)
colnames(System.Immune.Hydrolase.Sequences.df)<-c("PDB ID","Title","Residue Count","Sequence")

#---------------------------------------------------Sequence Alignments-------------------------------------------
#1
System.Immune.Hydrolase.Sequences.Alignment<- msa(System.Immune.Hydrolase.Sequences.df[,4],type="protein")
print(System.Immune.Hydrolase.Sequences.Alignment, show="complete")
#2
System.Immune.Hydrolase.Sequences.Alignment.Clustal <- msa(System.Immune.Hydrolase.Sequences.df[,4],                                                              type="protein", "ClustalW")
print(System.Immune.Hydrolase.Sequences.Alignment.Clustal, showConsensus=FALSE, halfNrow=3)
#3
System.Immune.Hydrolase.Sequences.Alignment.General<-msa(System.Immune.Hydrolase.Sequences.df[,4],
                                               type="protein",order="aligned")

System.Immune.Hydrolase.Sequences.Alignment.Convert<- msaConvert(System.Immune.Hydrolase.Sequences.Alignment, type="seqinr::alignment")
d <- dist.alignment(System.Immune.Hydrolase.Sequences.Alignment.Convert, "identity") 
Matrix.distance<-as.matrix(d)
dimnames(Matrix.distance)<-list(System.Immune.Hydrolase.Sequences.df[,2],System.Immune.Hydrolase.Sequences.df[,2])
System.Immune.Hydrolase.Tree <- nj(Matrix.distance)

msaPrettyPrint(System.Immune.Sequences.Classification.Alignment.Clustal, 
               output="tex", 
               y=c(1,100), 
               subset=c(1:5), 
               showNames="none", showLogo="top", 
               alFile ="data\\System_Immune_Hydrolase_Alignment_1.fasta",
               file="data\\System_Immune_Hydrolase_Alignment_1.tex",
               logoColors="rasmol", 
               shadingMode="functional",
               shadingModeArg="structure",
               showLegend=FALSE, 
               furtherCode=c("\\defconsensus{.}{lower}{upper}", "\\showruler{1}{top}",
               askForOverwrite=FALSE))

#----------------------------------------------------Tables------------------------------------------------------------------

Table.1<-xtable(System.Immune.Hydrolase.Sequences.df)

#------------Figures to be added by Students in the Classroom------------------------------------------------------------------
Figure.1<-plot(System.Immune.Hydrolase.Tree, main="Phylogenetic Tree of Immune System Hydrolase Sequences",cex=0.5)

#---------------------------------------------------References----------------------------------------------------------------

Reference.1<-c("Thompson, J. D., Higgins, D. G., and Gibson, T. J. (1994)",
                "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, 
                  osition-specific gap penalties and weight matrix choice.",
                "Nucleic Acids Res. 22(22):4673-4680. DOI: 10.1093/nar/22.22.4673.")

Reference.2<-c("Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K.", 
              "Li, W., Lopez, R., McWilliam, H., Remmert, M., Soeding, J., Thompson, J. D., and Higgins, D. G. (2011)", 
              "Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.", 
              "Mol. Syst. Biol. 7:539. DOI: 10.1038/msb.2011.75.")

Reference.3<-c("Beitz, E. (2000) TeXshade: shading and labeling of multiple sequence alignments using LaTeX2e.",
                "Bioinformatics 16(2):135-139.",
                "DOI: 10.1093/bioinformatics/16.2.135.")
#---------------------------------------------------Function Library----------------------------------------------------------

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

