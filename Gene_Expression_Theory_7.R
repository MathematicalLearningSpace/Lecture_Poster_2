#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------------------R API -----------------------------------------------------------------
library(xtable);
#-------------------------------------------------Normal Mode Analysis---------------------------------------------------
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);library(phangorn);library(proteomics);library(dcGOR)
library(MDplot);library(UniProt.ws);library(circlize);library(BioPhysConnectoR);library(ape)
#------------------------------------------------Arabidopsis Thaliana----------------------------------------------------
library(myTAI);library(seqinr);library(org.At.tair.db);library(dplyr);library(BioSeqClass);library(Rcpi);library(protr)
library(seqinr);library(Biostrings);library(Peptides);library(titrationCurves)
#------------------------------------------------Ontology--------------------------------------------------------------
library(dcGOR);library(ontologyIndex);library(ggplot2);library(tidyr);library(PearsonDS);library(CHNOSZ);library(stringi);library(stringr)
detach("package:Biostrings", unload=TRUE)
#------------------------------------------------Data------------------------------------------------------------
data("AAdata");data(thermo);data(aa.table)
thermo$element

Thermo.protein.df<-as.data.frame(thermo$protein)
Thermo.protein.df<-Thermo.protein.df[with(Thermo.protein.df,order(Thermo.protein.df$Cys,decreasing=TRUE)),]
pf <- as.data.frame(protein.formula(aa.strand))
#------------------------------------------------AA Data-----------------------------------------------------------

AA.df<-as.data.frame(aa.table)
AAData.df<-as.data.frame(AAdata)
AA.list<-unique(AA.df$aa1)
AA.df.mass<-AA.df$mass
AA.df<-AA.df[with(AA.df,order(AA.df$mass,decreasing=TRUE)),]

#------Ratios Designs for the Classroom------------------------------------------------------

#--------------------------------------FASTA Files for Classroom---------------------------
AT.Chromosome.3<-read.fasta(file ='Arabidopsis_thaliana.TAIR10.31.dna.chromosome.3.fa.gz')
AT.Chromosome.3[[1]][20339504:20341103]
AA.AT3G54890.3<-seq2aa("AT3G54890",AT.Chromosome.3[[1]][20339961:20340922])
#----------------------------Formatted Student Notes of FASTA Files------------------------
AT.Protein.Mitochondrion.1<-readFASTA('ArabidopsisThalianaProteinMitochondrion.txt')
AT.Protein.Mitochondrion.1.117<-AT.Protein.Mitochondrion.1[[117]]
AT.Protein.Mitochondrion.1.117.Motif.RRRR<-stri_locate_all_regex(AT.Protein.Mitochondrion.1.117, "(?=RRRR)")
#AT.Protein.Mitochondrion.1.117<-s2c(AT.Protein.Mitochondrion.1.117)
#------------------------------------------Properties--------------------------------------
AT.Protein.Mitochondrion.1.117.AA<-seq2aa("117",AT.Protein.Mitochondrion.1.117)
AT.Protein.Mitochondrion.1.117.mass<-mass(as.chemical.formula(protein.formula(AT.Protein.Mitochondrion.1.117.AA)))
AT.Protein.Mitochondrion.1.117.entropy<-entropy(as.chemical.formula(protein.formula(AT.Protein.Mitochondrion.1.117.AA)))
AT.Protein.Mitochondrion.1.117.composition<-aasum(AT.Protein.Mitochondrion.1.117.AA)
AT.Protein.Mitochondrion.1.117.Protein.Formula<-as.data.frame(protein.formula(AT.Protein.Mitochondrion.1.117.AA))
AT.Protein.Mitochondrion.1.117.Hydrogen.Carbon.ratio<-Hydrogen.Carbon.ratio(AT.Protein.Mitochondrion.1.117.Protein.Formula)
AT.Protein.Mitochondrion.1.117.Oxygen.Carbon.ratio<-Oxygen.Carbon.ratio(AT.Protein.Mitochondrion.1.117.Protein.Formula)
AT.Protein.Mitochondrion.1.117.Nitrogen.Carbon.ratio<-Nitrogen.Carbon.ratio(AT.Protein.Mitochondrion.1.117.Protein.Formula)
AT.Protein.Mitochondrion.1.117.Sulfur.Carbon.ratio<-Sulfur.Carbon.ratio(AT.Protein.Mitochondrion.1.117.Protein.Formula)
AT.Protein.Mitochondrion.1.117.carbon.Oxidation.average<-carbon.Oxidation.average(AT.Protein.Mitochondrion.1.117.Protein.Formula)

#---------------------------Table Design----------------------
AT.Protein.Mitochondrion.1.117.Statistics.df<-data.frame()
AT.Protein.Mitochondrion.1.117.Statistics.df<-rbind(c(AT.Protein.Mitochondrion.1.117.mass))
colnames(AT.Protein.Mitochondrion.1.117.Statistics.df)<-c("Total Mass")

#-------------Transformations to be designed by Students with this Design Pattern-------------------------------------------------
Transformation.1<-function(X)
  {
  Z<-NULL;
  w.1<-1
  for(i in 1:length(X))
    {
      Z[i]<-w.1.*X[i]
    }
  output<-list()
  output$X<-X
  output$Z<-Z
  return(output)
  }
X<-runif(100,0,1)
test.Transformation.1<-Transformation.1(X)
test.Transformation.1
#------------------------------------Tables-----------------------------------------------------------
Table.1<-xtable(AT.Protein.Mitochondrion.1.117.Statistics.df)

#---------------------------Figures for Presentation in Classroom----------------------------------------------------------
Figure.1<-plot(AA.df$mass,type = 'l', 
               col = "black", lwd = 2, lty=1,
               xaxt='n',
               main = "Amino Acid Mass Decreasing Sequence")
grid()
Axis(side=1, at=1:length(AA.df$aa3), 
     labels=AA.df$aa3,cex.axis=0.35)
legend("topright", col = c("black"), 
       lty = 1:1, cex=0.75,
       legend = c("Mass"))

Figure.2<-plot(Thermo.protein.df$Cys,type = 'l', 
               col = "black", lwd = 2, lty=1,
               xaxt='n',
               main = "Amino Acid CYS Decreasing Sequence and TYR")
grid()
Axis(side=1, at=1:length(Thermo.protein.df$Cys), 
     labels=Thermo.protein.df$abbrv,cex.axis=0.35)
lines(Thermo.protein.df$Tyr)
legend("topright", col = c("black"), 
       lty = 1:1, cex=0.75,
       legend = c("CYS","TYR"))

#----------Function Library to be Modified by Students in the Classroom-------------------------------------------------
#average state of oxidation
carbon.Oxidation.average<-function(d)
{
  molecule.charge=0
  coa<-(-d$H+3*d$N+2*d$O+2*d$S+molecule.charge)/d$C
  return(coa)
}
Hydrogen.Carbon.ratio <- function(d) d$H/d$C
Oxygen.Carbon.ratio <- function(d) d$O/d$C
Nitrogen.Carbon.ratio <- function(d) d$N/d$C
Sulfur.Carbon.ratio <- function(d) d$S/d$C

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


