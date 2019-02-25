#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API 
library(rcellminer);library(rcellminerData);library(bioCancer);library(trena);library(ctSGE);library(anamiR) 
library(srnadiff);library(RmiR);library(CancerSubtypes);library(GRENITS);library(compEpiTools);library(trigger);library(CVE);library(qgraph)
library(sgnesR);library(rcdk);library(sqldf)

#--------------------------------------------Data----------------------------------------------
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()

selectedOncogenes <- c("ABL1", "ALK", "BRAF", "CCND1", "CCND3", "CCNE1", "CCNE2", 
                       "CDC25A", "EGFR", "ERBB2", "EZH2", "FOS", "FOXL2", "HRAS", 
                       "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "MET", "MOS", "MYC", 
                       "NRAS", "PDGFB", "PDGFRA", "PIK3CA", "PIK3CB", "PIK3CD", 
                       "PIK3CG", "PIM1", "PML", "RAF1", "REL", "RET", "SRC", "STK11", 
                       "TP63", "WNT10B", "WNT4", "WNT2B", "WNT9A", "WNT3", "WNT5A", 
                       "WNT5B", "WNT10A", "WNT11", "WNT2", "WNT1", "WNT7B", "WISP1", 
                       "WNT8B", "WNT7A", "WNT16", "WISP2", "WISP3", "FZD5", "FZD1")
#---------------------------------------Function Library------------------------------------

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
