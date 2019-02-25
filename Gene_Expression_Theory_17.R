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
f.1<-function(gene.1,gene.2,Visualization=FALSE)
{
  acceptablePlotTypes <- c("drug", "cop", "exp", "xai", "exo", "mut", "mir", "pro", "mda")
  copPrefix <- "cop"
  expPrefix <- "exp"
  xaiPrefix <- "xai"
  exoPrefix <- "exo"
  mutPrefix <- "mut"
  mirPrefix <- "mir"
  proPrefix <- "pro"
  #------Filters------
  copData <- molData[[copPrefix]]
  expData <- molData[[expPrefix]]
  xaiData <- molData[[xaiPrefix]]
  exoData <- molData[[exoPrefix]]
  mutData <- molData[[mutPrefix]]
  mirData <- molData[[mirPrefix]]
  proData <- molData[[proPrefix]]
  
  gene.1<-gene.1
  gene.2<-gene.2
  gene<-c(gene.1,gene.2)
  
  # Get the cell lines names for cell lines meeting particular thresholds
  copKnockdown <- names(which(molData[["cop"]][paste0("cop", gene), ] < -1))
  expKnockdown <- names(which(molData[["exp"]][paste0("exp", gene), ] < -1.5))
  mutKnockdown <- names(which(molData[["mut"]][paste0("mut", gene), ] == 1))
  
  # Make composite pattern
  pattern <- rep(0, length(molData[["cop"]][paste0("cop", gene), ]))
  names(pattern) <- names(molData[["cop"]][paste0("cop", gene), ])
  tmp <- Reduce(union, list(copKnockdown, expKnockdown, mutKnockdown))
  pattern[tmp] <- 1
  
  #-----Summarize-----
  Table.1<-NULL
  
  #-----Visualize-----
  if(Visualization)
  {
    #------------------Profile Visualization------
    plots <- c("exp","exp") 
    plotCellMiner(drugAct, molData, plots, NULL, gene)
    # Composite plot data
    extraPlot <- list()
    extraPlot[["title"]] <- "Composite Pattern"
    extraPlot[["label"]] <- "Knockdown Composite (Binary)"
    extraPlot[["values"]] <- pattern
    
    plotCellMiner(molData=molData, plots=c("cop", "exp", "mut"), gene=gene, extraPlot=extraPlot)
    
  }
  
  #Object Model
  output<-list()
  output$Genes<-gene
  output$Table.1<-Table.1
  return(output)
  
  
}
test.f.1<-f.1("TP53","MDM2",Visualization=TRUE)
test.f.1

#-------------Function Template Library for Classroom Presentation and Modification---------------------
f.2<-function(X)
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
test.f.2<-f.2(letters)
test.f.2
