#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------------R API -----------------------------------------------------------
library(triplex);library(xtable);library(seqinr);library(Biostrings);library(BSgenome.Hsapiens.UCSC.hg19)
library(heatmaps);library(BSgenome);library(GenomeGraphs);library(biomaRt)
#--------------------------------Data--------------------------------------
strand.types<-c("Parallel-First-0","Parallel-First-1",
                "Parallel-Second-2","Parallel-Second-3",
                "AntiParallel-Second-4", "AntiParallel-Second-5",
                "AntiParallel-First-6", "AntiParallel-First-7")

genome <- BSgenome.Hsapiens.UCSC.hg19
chromosome.18<<-genome[["chr18"]]
k<-3
n<-(length(chromosome.18)/10^k)
t.18.3 <- triplex.search(chromosome.18[1:n],min_score=17,min_len=8)
t.18.3.Parallel.Second <- triplex.search(chromosome.18[1:n],
                                         type=c(2,3), min_score=10, p_value=1)
t.18.3.decreasing<-t.18.3[order(score(t.18.3), decreasing=TRUE)]

#-------------------------------Tables to be added by Students--------------------------------------
Table.1<-""
Table.2<-""
Table.3<-""

#-----Figures to be presented in the Classroom---------------------------------
par(mfrow = c(2,2))
Figure.1<-hist(t.18.3@elementMetadata@listData$type)
Figure.2<-hist(t.18.3@elementMetadata@listData$score)
Figure.3<-hist(t.18.3@ranges@width)
Figure.4<-hist(t.18.3@params)

Figure.5<-triplex.diagram(t.18.3[1])

Figure.6<-triplex.3D(t.18.3[1], A.col = "red", T.col = "brown", 
                     G.col = "green", C.col = "blue", bbone.col = "violet", 
                     bgr.col = "white", bbone.n = 20)

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

ideog <- makeIdeogram(chromosome = "18") 
genomeAxis <- makeGenomeAxis() 
genesplus <- makeGeneRegion(start = 0, end =n, strand = "+", chromosome = "18", biomart = mart) 
genesminus <- makeGeneRegion(start = 0, end =n, strand = "-", chromosome = "18", biomart = mart) 
tall <- as(t.18.3,"GRanges")
ta <- makeAnnotationTrack( start = start(tall[which(end(tall)<n)]), 
                           end = end(tall[which(end(tall)<n)]), 
                           feature = "gene_model", 
                           dp = DisplayPars(gene_model = "grey") 
                           ) 
ta1 <- makeAnnotationTrack( start = start(tall[which(end(tall)<n & score(tall)>17)]), 
                            end = end(tall[which(end(tall)<n & score(tall)>17)]), 
                            feature = "gene_model", 
                            dp = DisplayPars(gene_model = "darkblue")
                            ) 
hr1 <- makeRectangleOverlay(start = 11000, end = 13000)
hr2 <- makeRectangleOverlay(start = 45000, end = 50000)
roR <- makeRectangleOverlay(start = .1, end = .3, 
                            coords = "absolute",
                            dp = DisplayPars(fill = "grey", alpha = .2, lty = "dashed"),
                            region = c(.4,.7))
Figure.7<-gdPlot(list(fwd = genesplus, 
            GENES = genomeAxis, 
            rev = genesminus, 
            weak = ta, strong = ta1, chr18 = ideog), 
            minBase = 0, 
            maxBase =n, labelRot = 0,overlays = list(hr1, hr2,roR))

#-------------References to be added by Students in the Classroom -----------------------------

Reference.1<-c("",
               "",
               "")

#-----Function Library to be added by Students in the Classroom--------------------------

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
