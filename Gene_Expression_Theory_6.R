#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(stringr);library(stringi);library(chinese.misc)

#-------------------------------------------Data----------------------------------------------------------
#---------------HGNC Gene Family Data Set For Classroom --------------------------------------------------
Gene.Family.HGNC <- read_delim("HGNC_Gene_Family_Data.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
View(Gene.Family.HGNC)
Gene.Family.HGNC$Chromosome
#------------------------------------------Choose a Chromosome --------------------------------------------
#--------------------------16 Percent Cancer and 8p + 15 high mutation rate (Point and Missense)-------------
#------------------------------------------Wikipedia--------------------------------------------------------
pattern<-"8q.*"
Gene.Chromosome.8.q<-match_pattern(pattern, Gene.Family.HGNC$Chromosome)
Gene.Chromosome.8.q<-Gene.Family.HGNC$`Approved Symbol`[grep('8q.*',Gene.Family.HGNC$Chromosome)]

pattern.2<-"WDR*"
Gene.WDR.Chromosomes<-match_pattern(pattern.2, Gene.Family.HGNC$`Approved Symbol`)
Gene.WDR.Chromosomes<-Gene.Family.HGNC$Chromosome[grep('WDR*',Gene.Family.HGNC$`Approved Symbol`)]

pattern.3<-"8p22.*"
Gene.Chromosome.8.p.22<-match_pattern(pattern.3, Gene.Family.HGNC$Chromosome)
Gene.Chromosome.8.p.22<-Gene.Family.HGNC$`Approved Symbol`[grep('8p22.*',Gene.Family.HGNC$Chromosome)]

pattern.4<-"8p23.*"
Gene.Chromosome.8.p.23<-match_pattern(pattern.4, Gene.Family.HGNC$Chromosome)
Gene.Chromosome.8.p.23<-Gene.Family.HGNC$`Approved Symbol`[grep('8p23.*',Gene.Family.HGNC$Chromosome)]

Gene.Chromosome.8.p.23.study<-c("DEFA*","OR7E*","ZNF*","MIR*")

Transformations<-c("translocations","inversions","deletions","duplications")
Defensins.publications <- read_csv("Defensins/publications.csv")
View(publications)
#-----------------------------------------Regular Expressions---------------------------------------------

Gene.DEFA<-grep("DEFA",Gene.Chromosome.8.p.23)
Gene.OR7E<-grep("OR7",Gene.Chromosome.8.p.23)
Gene.ZNF<-grep("ZNF",Gene.Chromosome.8.p.23)
Gene.MIR<-grep("MIR",Gene.Chromosome.8.p.23)

#-----------------------------------------Peptide Transformations------------------------------------------------
Transformations<-c("translocations","inversions","deletions","duplications")


Chromosome_Analysis.df<-data.frame()
Chromosome_Analysis.df<-rbind(TF_1(Chromosome.Study))
colnames(Chromosome_Analysis.df)<-c("Transformation")
rownames(Chromosome_Analysis.df)<-c("Chromosome_1")

#------------------------------------------Definitions----------------------------------------------------
#----------------------------------Formatted Student Notes------------------------------------------------
Definitions.df<-as.data.frame(read_delim("Gene_Expression_Definitions.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE))

#------------------------------------------Theorems-------------------------------------------------------
#----------------------------------Formatted Student Notes------------------------------------------------
Theorems.df<-as.data.frame(read_delim("Gene_Expression_Theorems.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE))


#------------------------------------------Tables---------------------------------------------------------
Table.1<-xtable(Chromosome_Analysis.df)
Table.2<-xtable(Definitions.df)
Table.3<-xtable(Theorems.df)

#------------------------------------------Figures---------------------------------------------------------


#------------------------------------------Function Library-----------------------------------------------

TF_1->function(x)
{
  TF.data.local<-NULL
  #Design of Transformations goes here
  return(list(TF.Local=TF.data.local))
  
}
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
