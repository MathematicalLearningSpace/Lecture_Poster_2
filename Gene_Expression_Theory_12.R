#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------R API ---------------------------------------
library(xtable);library(RJaCGH);library(GLAD);library(seqCNA)
#----------------------------------Data for Classroom-----------------------------------------
data(snijders) ;data(cytoband);data(arrayCGH)

array <- list(arrayValues=array2, 
              arrayDesign=c(4,4,21,22))
class(array) <- "arrayCGH"

profileCGH <- as.profileCGH(gm07408)
res <- glad(profileCGH, mediancenter=FALSE,
            smoothfunc="lawsglad", bandwidth=10, round=1.5,
            model="Gaussian", lkern="Exponential", qlambda=0.999,
            base=FALSE,
            lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
            type="tricubic", param=c(d=6),
            alpha=0.001, msize=2,
            method="centroid", nmax=8,
            verbose=FALSE)

res.df<-as.data.frame(res)
res$SigmaC
res$BkpInfo

gm07408.y <- gm07408$LogRatio[!is.na(gm07408$LogRatio)] 
gm07408.Pos <- gm07408$PosBase[!is.na(gm07408$LogRatio)] 
gm07408.Chrom <- gm07408$Chromosome[!is.na(gm07408$LogRatio)]

id <- order(gm07408.Chrom, gm07408.Pos) 
y <- y[id] 
Pos <- Pos[id] 
Chrom <- Chrom[id]

#----------------------------------MC Models-----------------------------------

jump.parameters <- list(sigma.tau.mu = rep(0.01, 4), 
                        sigma.tau.sigma.2 = rep(0.05, 4), 
                        sigma.tau.beta = rep(0.1, 4), 
                        tau.split.mu = 0.1)

fit <- RJaCGH(y = y, Pos = Pos, Chrom = Chrom, 
              model = "genome", var.equal = TRUE, k.max = 4, 
              burnin = 100, TOT = 100, 
              jump.parameters = jump.parameters, NC = 2, deltaT = 0.5)

fit.2 <- relabelStates(fit, window = 1.25)

names(fit)
length(fit[["array1"]]$fit.k)

summary.HMM.1 <- summary(fit) 
summary.HMM.1

summary(fit.2)

#---------------------------------Chromosome Models---------------------------

y2 <- gm07408$LogRatio[!is.na(gm07408$LogRatio)] 
Pos2 <- gm07408$PosBase[!is.na(gm07408$LogRatio)] 
Chrom2 <- gm07408$Chromosome[!is.na(gm07408$LogRatio)] 

id <- order(Chrom2, Pos2) 
y2 <- y2[id] 
Pos2 <- Pos2[id] 
Chrom2 <- Chrom2[id] 

fit.chrom <- RJaCGH(y = y2, Pos = Pos2, 
                    Chrom = Chrom2, 
                    model = "Chrom", k.max = 4, burnin = 100, TOT = 100, NC = 2, deltaT = 0.5)

summary(fit.chrom, array = "array1", Chrom = 8)

sequence <- states(fit.chrom) 
sequence.averaged <- modelAveraging
sequence.averaged <- modelAveraging(fit.chrom)
s.means <- smoothMeans(fit.chrom)

#-------------------------------Common Regions--------------------------------
gm01750LR <- gm01750$LogRatio 
gm04435LR <- gm04435$LogRatio 
gm05296LR <- gm05296$LogRatio 
gm07408LR <- gm07408$LogRatio 

not.NA <- !is.na(gm01750LR) & !is.na(gm04435LR) & !is.na(gm05296LR) & !is.na(gm07408LR) 
gm01750LR <- gm01750LR[not.NA] 
gm04435LR <- gm04435LR[not.NA] 
gm05296LR <- gm05296LR[not.NA]
gm07408LR <- gm07408LR[not.NA] 

Pos8 <- gm07408$PosBase[not.NA] 
Chrom8 <- gm07408$Chromosome[not.NA]

id <- order(Chrom8, Pos8) 
Pos8 <- Pos8[id] 
Chrom8 <- Chrom8[id] 

fit.arrays <- RJaCGH(y = cbind(gm01750LR = gm01750LR[id], 
                               gm04435LR = gm04435LR[id],
                               gm05296LR = gm05296LR[id],
                               gm07408LR = gm07408LR[id]), Pos = Pos8, 
                     Chrom = Chrom8, model = "genome", k.max = 4, burnin = 100, 
                     TOT = 100, NC = 2, deltaT = 0.5)
summary(fit.arrays)
summary(fit.arrays, array = "gm07408LR")

Regions.Gain <- pREC_A(fit.arrays, p = 0.33, alteration = "Gain") 
Regions.Loss <- pREC_A(fit.arrays, p = 0.33, alteration = "Loss")
Regions.Gain

RG.df <- as.data.frame(print(Regions.Gain))
Regions <- pREC_S(fit.arrays, p = 0.75, alteration = "Gain", freq.array = 2)

#---------------------------------Diagnostics---------------------------------

maxK <- as.numeric(names(which.max(table(fit.2[["array1"]]$k)))) 

fit.2[["array1"]]$fit.k[[maxK]]$prob.mu
fit.2[["array1"]]$fit.k[[maxK]]$prob.sigma.2
fit.2[["array1"]]$fit.k[[maxK]]$prob.beta

fit.2[["array1"]]$prob.b
fit.2[["array1"]]$prob.d
fit.2[["array1"]]$prob.s
fit.2[["array1"]]$prob.c

#---------------------------------Tables---------------------------------------

Table.1<-xtable(res.df[grep("^8",res.df$Chromosome),])

round(fit[["array1"]]$fit.k[[8]]$state.labels, 3)
seq.states <- modelAveraging(fit.arrays, array = "gm07408LR")[["gm07408LR"]]$states 

Table.2<-xtable(table(seq.states, gm07408$Statut[not.NA]))
Table.3<-xtable(RG.df)

#----------Figures to be presented in the Classroom--------------------------

Figure.1<-plotProfile(res, unit=3, Bkp=TRUE, labels=FALSE,
            Smoothing="Smoothing", plotband=TRUE, cytoband = cytoband)

text <- list(x=c(90000,200000),y=c(0.15,0.3),labels=c("NORMAL","GAIN"), cex=2)
Figure.2<-plotProfile(res, unit=3, Bkp=TRUE, labels=TRUE, Chromosome=8,
            Smoothing="Smoothing", text=text, main="Chromosome 8", cytoband = cytoband)

Figure.3<-arrayPlot(array,"Log2Rat", main="Spatial Image of CGH")

Figure.4<-plot(fit, array = "array1", cex = 1.1)
Figure.5<-plot(relabelStates(fit, window = 0.25)) 
Figure.6<-plot(relabelStates(fit, window = 2))
Figure.7<-trace.plot(fit.2, array = "array1")

Figure.8<-plot(fit.chrom, array="array1")
Figure.9<-plot(fit.chrom,array="array1",Chrom=8)
Figure.10<-genomePlot(fit.chrom)

Figure.11<-plot(fit.arrays, show = "frequency")
Figure.12<-plot(fit.arrays, show = "average")
Figure.13<-plot(Regions, cex.axis = 0.6)

Figure.14<-hist(y)
Figure.15<-hist(Chrom8)

#-------------------------------References------------------------------------
Reference.1<-c("Rueda OM, Diaz-Uriarte R", 
  "Flexible and accurate detection of genomic copy-number changes from aCGH.", 
  "PLoS Comput Biol. 2007, 3(6):1115-1122.") 

Reference.2<-c("McCarroll SA, Altshuler DM", 
  "Copy-number variation and association studies of human disease.", 
  "Nat Genet 2007, 39(7 Suppl):S37-S42, [http://dx.doi.org/10.1038/ng2080].")


#-------------------------------Function Library-----------------------------
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
