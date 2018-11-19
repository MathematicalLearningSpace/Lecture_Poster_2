#--------------------------------------R API -----------------------------------------
library(readr);library(xtable)

#---------------------------------------Data from Pharos (https://pharos.nih.gov)---------------------------------------------
#---------------------------------------Student build the csv files in the classroom-----------------------------------------
#---------------------------------------Cancer Research Publications---------------------------------------------------------

AML_Cancer.publications <- read_csv("AML/publications.csv")
#View(AML_Cancer.publications)

Bladder_Cancer.publications <- read_csv("BladderCancer/publications.csv")
#View(Bladder_Cancer.publications)

Breast_Cancer.publications <- read_csv("Breast Cancer/publications.csv")
#View(Breast_Cancer.publications)

CML_Cancer.publications <- read_csv("CML/publications.csv")
#View(CML_Cancer.publications)

Colorectal_Cancer.publications <- read_csv("Colorectal/publications.csv")
#View(Colorectal_Cancer.publications)

Endometrial_Cancer.publications <- read_csv("Endometrial/publications.csv")
#View(Endometrial_Cancer.publications)

Gastric_Cancer.publications <- read_csv("GastricCancer/publications.csv")
#View(Gastric_Cancer.publications)

Hepatocellular_Cancer.publications <- read_csv("Hepatocellular/publications.csv")
#View(Hepatocellular_Cancer.publications)

NonsmallCellLung_Cancer.publications <- read_csv("NonsmallCellLung/publications.csv")
#View(NonsmallCellLung.publications)

Pancreatic_Cancer.publications <- read_csv("Pancreatic/publications.csv")
#View(Pancreatic_Cancer.publications)

Prostate_Cancer.publications <- read_csv("Prostate/publications.csv")
#View(Prostate_Cancer.publications)

smallCellLung_Cancer.publications <- read_csv("smallCellLung/publications.csv")
#View(smallCellLung_Cancer.publications)

Thyroid_Cancer.publications <- read_csv("Thyroid Cancer/publications.csv")
#View(Thyroid_Cancer.publications)

Cancer.publications<-c(AML_Cancer.publications,Bladder_Cancer.publications, Breast_Cancer.publications,
                       CML_Cancer.publications, Colorectal_Cancer.publications, Endometrial_Cancer.publications,
                       Gastric_Cancer.publications, Hepatocellular_Cancer.publications, NonsmallCellLung_Cancer.publications,
                       Pancreatic_Cancer.publications, Prostate_Cancer.publications, smallCellLung_Cancer.publications,
                       Thyroid_Cancer.publications )

#-------------------------------------Biological Signaling Networks---------------------------------------------

Apoptosis_Signaling.publications <- read_csv("Apoptosis/publications.csv")
View(Apoptosis_Signaling.publications)

AR_Signaling.publications <- read_csv("AR/publications.csv")
View(AR_Signaling.publications)

Calcium_Signaling.publications <- read_csv("Calcium/publications.csv")
View(Calcium_Signaling.publications)

Cell_Cycle_G1_S_Signaling.publications <- read_csv("Cell_Cycle_G1/publications.csv")
View(Cell_Cycle_G1_S_Signaling.publications)

HH_Signaling.publications <- read_csv("HH/publications.csv")
View(HH_Signaling.publications)

HIF_1_Signaling.publications <- read_csv("HIF_1/publications.csv")
View(HIF_1_Signaling.publications)

JAK_STAT_Signaling.publications <- read_csv("JAK_STAT/publications.csv")
View(JAK_STAT_Signaling.publications)

KEAP1_NRF2_Signaling.publications <- read_csv("KEAP1_NRF2/publications.csv")
View(KEAP1_NRF2_Signaling.publications)

MAPK_ERK_Signaling.publications <- read_csv("MAPK_ERK/publications.csv")
View(MAPK_ERK_Signaling.publications)

NOTCH_Signaling.publications <- read_csv("NOTCH/publications.csv")
View(NOTCH_Signaling.publications)

Other_RAS_Signaling.publications <- read_csv("Other_RAS/publications.csv")
View(Other_RAS_Signaling.publications)

PI3K_Signaling.publications <- read_csv("PI3K/publications.csv")
View(PI3K_Signaling.publications)

Telomerase_Signaling.publications <- read_csv("Telomerase/publications.csv")
View(Telomerase_Signaling.publications)

TGFB_Signaling.publications <- read_csv("TGFB/publications.csv")
View(TGFB_Signaling.publications)

WNT_Signaling.publications <- read_csv("WNT/publications.csv")
View(WNT_Signaling.publications)

Network.Signaling.publications<-c(Apoptosis_Signaling.publications,AR_Signaling.publications,Calcium_Signaling.publications,
                                  Cell_Cycle_G1_S_Signaling.publications,HH_Signaling.publications,HIF_1_Signaling.publications,
                                  JAK_STAT_Signaling.publications,KEAP1_NRF2_Signaling.publications,MAPK_ERK_Signaling.publications, 
          NOTCH_Signaling.publications, Other_RAS_Signaling.publications, PI3K_Signaling.publications, 
          Telomerase_Signaling.publications,TGFB_Signaling.publications, WNT_Signaling.publications)


#-----------------------------------------------------------------------------------------------------------

Network.signaling.df<-data.frame()

Network.signaling.df<-cbind(length(Apoptosis_Signaling.publications),
length(AR_Signaling.publications),
length(Calcium_Signaling.publications),
length(Cell_Cycle_G1_S_Signaling.publications),
length(HH_Signaling.publications),
length(HIF_1_Signaling.publications),
length(JAK_STAT_Signaling.publications),
length(KEAP1_NRF2_Signaling.publications),
length(MAPK_ERK_Signaling.publications), 
length(NOTCH_Signaling.publications),
length(Other_RAS_Signaling.publications),
length(PI3K_Signaling.publications),
length(Telomerase_Signaling.publications),
length(TGFB_Signaling.publications),
length(WNT_Signaling.publications))
colNames(Network.signaling.df)<-c(Apoptosis_Signaling,
                                  AR_Signaling,
                                  Calcium_Signaling,
                                  Cell_Cycle_G1_S_Signaling,
                                  HH_Signaling,
                                  HIF_1_Signaling,
                                  JAK_STAT_Signaling,
                                  KEAP1_NRF2_Signaling,
                                  MAPK_ERK_Signaling,
                                  NOTCH_Signaling,
                                  Other_RAS_Signaling,
                                  PI3K_Signaling,
                                  Telomerase_Signaling,
                                  TGFB_Signaling,
                                  WNT_Signaling)


Cancer.publications.df<-data.frame()
Cancer.publications.df<-cbind(length(AML_Cancer.publications$`Pubmed ID`),
                              length(Bladder_Cancer.publications$`Pubmed ID`), 
                              length(Breast_Cancer.publications$`Pubmed ID`),
                              length(CML_Cancer.publications$`Pubmed ID`), 
                              length(Colorectal_Cancer.publications$`Pubmed ID`), 
                              length(Endometrial_Cancer.publications$`Pubmed ID`),
                              length(Gastric_Cancer.publications$`Pubmed ID`), 
                              length(Hepatocellular_Cancer.publications$`Pubmed ID`), 
                              length(NonsmallCellLung_Cancer.publications$`Pubmed ID`),
                              length(Pancreatic_Cancer.publications$`Pubmed ID`), 
                              length(Prostate_Cancer.publications$`Pubmed ID`), 
                              length(smallCellLung_Cancer.publications$`Pubmed ID`),
                              length(Thyroid_Cancer.publications$`Pubmed ID` ))

colNames(Cancer.publications.df)<-c("AML_Cancer",
                                    "Bladder_Cancer",
                                    "Breast_Cancer",
                                    "CML_Cancer", 
                                    "Colorectal_Cancer",
                                    "Endometrial_Cancer",
                                    "Gastric_Cancer",
                                    "Hepatocellular_Cancer",
                                    "NonsmallCellLung_Cancer",
                                    "Pancreatic_Cancer",
                                    "Prostate_Cancer",
                                    "smallCellLung_Cancer",
                                    "Thyroid_Cancer")

#----------------------------------------------Tables------------------------------------------------------

Table.1<-xtable(Cancer.publications.df)
Table.2<-xtable(Network.signaling.df)

#----------------------Figures to be added by Students-----------------------------------------------------



#---------------------------------------------Function Library---------------------------------------------

