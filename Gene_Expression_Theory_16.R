#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(text2vec);library(data.table);library(magrittr);library(glmnet);library(rentrez);library(readr);library(LDAvis)
library(GSVAdata);library(org.Hs.eg.db);library(GSEABase);library(GSAR);library(MASS)
#------------------------------------------Data-----------------------------------------------

data(p53DataSet) #p53 Dataset of the NCI-60 Cell Lines
data(c2BroadSets) #C2 collection of canonical pathways from MSigDB 3.0
data(leukemia) #Leukemia Data by Armstrong et al. (2002) from the Broad Institute

C2 <- as.list(geneIds(c2BroadSets))
len <- length(C2)

genes.entrez <- unique(unlist(C2))
genes.symbol <- array("",c(length(genes.entrez),1)) 
symbol.map <- org.Hs.egSYMBOL 
mapped.genes <- mappedkeys(symbol.map)
symbol.map.list <- as.list(symbol.map[mapped.genes]) 

for (ind in 1:length(genes.entrez)){ 
  if (length(symbol.map.list[[genes.entrez[ind]]])!=0) 
  genes.symbol[ind] <- symbol.map.list[[genes.entrez[ind]]] 
}

genes.no.mapping <- which(genes.symbol == "") 
if(length(genes.no.mapping) > 0){ 
  genes.entrez <- genes.entrez[-genes.no.mapping] 
  genes.symbol <- genes.symbol[-genes.no.mapping] 
}
#---------------------------------------Gene Lists---------------------------------------------
names(genes.symbol) <- genes.entrez 
p53genes <- rownames(p53DataSet) 
remained <- array(0,c(1,len))
for (k in seq(1, len, by=1)) { 
  remained[k] <- sum((genes.symbol[C2[[k]]] %in% p53genes) & 
                       (C2[[k]] %in% genes.entrez)) 
}
#---------------------------------------Pathways----------------------------------------------

C2 <- C2[(remained>=10)&&(remained<=500)] 
pathway.names <- names(C2)
c2.pathways <- list() 
for (k in seq(1, length(C2), by=1)) { 
  selected.genes <- which(p53genes %in% genes.symbol[C2[[k]]]) 
c2.pathways[[length(c2.pathways)+1]] <- p53genes[selected.genes] 
} 
names(c2.pathways) <- pathway.names
path.index <- which(names(c2.pathways) == "LU_TUMOR_VASCULATURE_UP")
target.pathway <- p53DataSet[c2.pathways[["LU_TUMOR_VASCULATURE_UP"]],] 
group.label <- c(rep(1,17), rep(2,33)) 

#----------------------------------Testing---------------------------------------------------
#---Multivariate Generalization of the Wald-Wolfowitz Runs Test
WW.pvalue <- WWtest(target.pathway, group.label) 
#---Multivariate Kolmogorov-Smirnov Test of Means
KS.pvalue <- KStest(target.pathway, group.label) 
#---Multivariate Mean Deviation Test of Means
MD.pvalue <- MDtest(target.pathway, group.label) 
#---Multivariate Radial Kolmogorov-Smirnov Test of Variance
RKS.pvalue <- RKStest(target.pathway, group.label) 
#---Multivariate Radial Mean Deviation Test of Variance
RMD.pvalue <- RMDtest(target.pathway, group.label) 
F.pvalue <- AggrFtest(target.pathway, group.label) 
GSNCA.pvalue <- GSNCAtest(target.pathway, group.label)

path.index.2 <- which(names(c2.pathways) == "LU_TUMOR_VASCULATURE_DN")
target.pathway.2 <- p53DataSet[c2.pathways[["LU_TUMOR_VASCULATURE_DN"]],]

WW.2.pvalue <- WWtest(target.pathway.2, group.label) 
KS.2.pvalue <- KStest(target.pathway.2, group.label) 
MD.2.pvalue <- MDtest(target.pathway.2, group.label) 
RKS.2.pvalue <- RKStest(target.pathway.2, group.label) 
RMD.2.pvalue <- RMDtest(target.pathway.2, group.label) 
F.2.pvalue <- AggrFtest(target.pathway.2, group.label) 
GSNCA.2.pvalue <- GSNCAtest(target.pathway.2, group.label) 

GSEA.df<-data.frame()
GSEA.df<-rbind(c("UP",WW.pvalue, KS.pvalue,MD.pvalue,RKS.pvalue,RMD.pvalue,F.pvalue,GSNCA.pvalue),
               c("Down",WW.2.pvalue, KS.2.pvalue,MD.2.pvalue,RKS.2.pvalue,RMD.2.pvalue,F.2.pvalue,GSNCA.2.pvalue))
colnames(GSEA.df)<-c("Direction","WW","KS","MD","RKS","RMD","F","GSNCA")

#----------------Gene Sets Net Correlations Analysis------------------------------------------

ngenes <- 20
nsamples <- 40
zero.vector <- array(0,c(1,ngenes))
cov.mtrx1 <- diag(ngenes)
cov.mtrx1[!diag(ngenes)] <- 0.1
cov.mtrx2 <- diag(ngenes)
cov.mtrx2[!diag(ngenes)] <- 0.1
mask <- diag(ngenes/4)
mask[!diag(ngenes/4)] <- 0.6
cov.mtrx2[1:(ngenes/4),1:(ngenes/4)] <- mask

gp1 <- mvrnorm((nsamples/2), zero.vector, cov.mtrx1)
gp2 <- mvrnorm((nsamples/2), zero.vector, cov.mtrx2)
gp <- rbind(gp1,gp2)
dataset <- aperm(gp, c(2,1))
rownames(dataset) <- as.character(c(1:ngenes))
sample.labels <- c(rep(1,20),rep(2,20))

pvalue <- GSNCAtest(object=dataset, group=sample.labels) 

geneSets <- list("set1"=as.character(c(1:20)), 
                 "set2"=as.character(c(11:40)), 
                 "set3"=as.character(c(31:40)))

results.1 <- TestGeneSets(object=dataset, group=sample.labels, 
                        geneSets=geneSets, test="KStest")

results.2 <- TestGeneSets(object=p53DataSet, 
                        group=group.label, 
                        geneSets=c2.pathways[1:3], 
                        min.size=10, 
                        max.size=100, test="GSNCAtest") 

#---Aggregated F-Test of Variance Using Fisher's Probability Combining Method
pvalue <- AggrFtest(object=dataset, group=sample.labels)

#-------------------------Minimum Spanning Trees---------------------------------------
dataset.2<-p53DataSet[1:40,1:20]
#--------High Directed Preorder Ranking of MST----------------------------------------
Wmat <- as.matrix(dist(dataset.2, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
#--------create a weighted undirectional graph from the weight matrix------------------
gr <- graph_from_adjacency_matrix(Wmat, weighted = TRUE, mode = "undirected")
#--------find the minimum spanning tree------------------------------------------------
MST <- mst(gr)
HDP.ranks <- HDP.ranking(MST)
radial.ranks <- radial.ranking(MST)
#---------Union of the First and Second Minimum Spanning Trees-------------------------
trees <- findMST2(dataset.2)
#---------Union of the First and Second Minimum Spanning Trees for PPI Networks-------
mst2.ppi <- findMST2.PPI(object=trees, return.MST2only=TRUE)
degs <- degree(mst2.ppi)
ind <- which(degs > 10)
V(mst2.ppi)$color <- "green"
V(mst2.ppi)$color[ind] <- "red"
#---------------------------------Gene Lists-----------------------------------

Gene.Set.UP.df<-data.frame()
Gene.Set.UP.df<-cbind(c2.pathways[1:3]$NAKAMURA_CANCER_MICROENVIRONMENT_UP)
colnames(Gene.Set.UP.df)<-c("Up")

Gene.Set.Down.df<-data.frame()
Gene.Set.Down.df<-cbind(c2.pathways[1:3]$NAKAMURA_CANCER_MICROENVIRONMENT_DN)
colnames(Gene.Set.Down.df)<-c("Down")


#---------------------------------Tables--------------------------------------------------------

Table.1<-xtable(GSEA.df)
Table.1
Table.2<-xtable(Gene.Set.UP.df)
Table.2
Table.3<-xtable(Gene.Set.Down.df)
Table.3

#---------------------------------Figures--------------------------------------------------------

par(mfcol = c(1, 1))
Figure.1<-plot(target.pathway[,1], xaxt = "n",type="l")
axis(1, 1:length(target.pathway[,1]),names(target.pathway[,1]),col.axis="purple",cex.axis=0.5)
#-----------------Add Additional Figures----------------------------
for(i in 2:ncol(target.pathway))
{
  if(i > 17)
    lines(target.pathway[,i],lty=i,col="red")
  else{
    lines(target.pathway[,i],lty=i,col="green")
  }
}
par(mfcol = c(1, 1))
Figure.2<-plot(target.pathway.2[,1], xaxt = "n",type="l",ylim=c(2,8))
axis(1, 1:length(target.pathway.2[,1]),names(target.pathway.2[,1]),col.axis="purple",cex.axis=0.5)
for(i in 2:ncol(target.pathway.2))
{
  #-----------------Add Additional Figures----------------------------
  if(i > 17)
    lines(target.pathway.2[,i],lty=i,col="red")
  else{
    lines(target.pathway.2[,i],lty=i,col="green")
  }
}

Figure.3<-plotMST2.pathway(object=p53DataSet[c2.pathways[[path.index]],], 
                 group=c(rep(1,17), rep(2,33)), name="LU_TUMOR_VASCULATURE_UP", 
                 legend.size=1.2, leg.x=-1.2, leg.y=2, 
                 label.size=1, label.dist=0.8, cor.method="pearson")

Figure.4<-plotMST2.pathway(object=p53DataSet[c2.pathways[[path.index.2]],], 
                           group=c(rep(1,17), rep(2,33)), name="LU_TUMOR_VASCULATURE_DN", 
                           legend.size=1.2, leg.x=-1.2, leg.y=2, 
                           label.size=1, label.dist=0.8, cor.method="pearson")

#-------------------------Figure Group 1-------------------------
par(mfcol = c(2, 2))
Figure.5<-plot(gr)
Figure.6<-plot(MST)
Figure.7<-plot(trees)
Figure.8<-plot(mst2.ppi)
#---------------------------------References----------------------------------------------------

Reference.1<-c("Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012)", 
                "Gene set analysis for self-contained tests: complex null and specific alternative hypotheses.", 
                "Bioinformatics 28, 3073-3080")
Reference.2<-c("Gene Set Enrichment Analysis",
               "http://software.broadinstitute.org/gsea/datasets.jsp")

#---------------------------------Function Library-----------------------------------------------
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
