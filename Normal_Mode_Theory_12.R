library(xtable)
library(Rknots)
library(bio3d)
library(igraph)
library(Rpdb)
library(DNAtools)
library(stringr)
library(stringi)
library(sfsmisc)
library(readr)
library(graphics)

#------------------Data------------------------------------------------

Results_RCSB_DNA_Ligases_HomoSapiens <- read_csv("Results_RCSB_ DNA_Ligases_HomoSapiens.csv")
#View(Results_RCSB_DNA_Ligases_HomoSapiens)

pdb<-read.pdb("1x9n.pdb")
pdb.2<-read.pdb("1fu1.pdb")

#-------------Classify Surface based on Temperature Factor--------------

test.classify.atoms<-classify.atoms(pdb$atoms$temp)
test.classify.atoms.2<-classify.atoms(pdb.2$atoms$temp)

#-------------Fit and predict with a cubic smooth spline-----------------------------------

cs <- smooth.spline(seq(1:665),grep("S",test.classify.atoms$Surface.loc))
cs
cs.predict<-predict(cs)

#-------------Numerical derivatives,i.e. First and Second---------------

d12<-D1D2(seq(1:665),grep("S",test.classify.atoms$Surface.loc))

#------------Approximate Integration------------------------------------

int.value<-integrate.xy(seq(1:665),grep("S",test.classify.atoms$Surface.loc))

#------------Test the Differences Between Factors--------------------------------------

t<-t.test(pdb$atoms$temp,pdb.2$atoms$temp, alternative = "g")
var.t<-var.test(pdb$atoms$temp,pdb.2$atoms$temp)
wilcox<-wilcox.test(pdb$atoms$temp, pdb.2$atoms$temp, alternative = "g")
ks<-ks.test(pdb$atoms$temp, pdb.2$atoms$temp, alternative = "l")

test.df<-data.frame()
test.df<-rbind(c(t$alternative, t$statistic,t$p.value),
               c(var.t$alternative, var.t$statistic,var.t$p.value),
               c(wilcox$alternative, wilcox$statistic,wilcox$p.value),
               c(ks$alternative, ks$statistic,ks$p.value))
colnames(test.df)<-c("Alternative","Statistic","p")
rownames(test.df)<-c("T","Variance", "Wilcox", "KS")

#------------------------Tables----------------------------------------

Table.1<-xtable(test.df)

table(test.classify.atoms$Surface.loc)
table(test.classify.atoms.2$Surface.loc)

#-------------------------Figures---------------------------------------
par(mfcol = c(2, 1))
Figure.1<-plot(pdb$atoms$temp,col="green")
Figure.1A<-plot(pdb.2$atoms$temp,col="blue")

par(mfcol = c(2, 2))
Figure.2<-hist(pdb$atoms$temp)
Figure.2A<-hist(pdb.2$atoms$temp)
Figure.2B<-plot(ecdf(pdb$atoms$temp), xlim = range(c(pdb$atoms$temp, pdb.2$atoms$temp)))
Figure.2C<-plot(ecdf(pdb.2$atoms$temp), add = FALSE, lty = "dashed")

par(mfcol = c(2, 1))
Figure.3<-plot(grep("S",test.classify.atoms$Surface.loc))
lines(cs <- smooth.spline(seq(1:665),grep("S",test.classify.atoms$Surface.loc)),col="blue")
lines(cs.predict<-predict(cs),col="red")
Figure.3A<-plot(grep("S",test.classify.atoms.2$Surface.loc))
lines(ss.2 <- smooth.spline(seq(1:2768),grep("S",test.classify.atoms.2$Surface.loc)),col="green")

Figure.4<-plot(d12$x,d12$D1)
lines(d12$x,d12$D2,col="Red")
lines(difference.quotient(seq(1:665),grep("S",test.classify.atoms$Surface.loc)),col="blue")

hc<-hclust(dist(pdb$atoms$temp))
dend1 <- as.dendrogram(hc)
dend2 <- cut(dend1, h = 0.3)
op <- par(mfrow =  c(2,2), mar = c(5,2,1,4))
plot(dend2$upper,horiz = FALSE,edgePar = list(col = "blue", lwd = 0.5),leaflab="none")
plot(dend2$upper, nodePar = list(pch = c(1,NA), cex = 0.8, lab.cex = 0.2),
     type = "t", center = TRUE,leaflab="none")
plot(dend2$upper, edgePar = list(col = 1:2, lty = 2:3,lab.cex = 0.2),
     dLeaf = 1, edge.root = TRUE,leaflab="none")
plot(dend2$upper, nodePar = list(pch = 2:1, cex = .4*2:1, col = 2:3, lab.cex = 0.2),
     horiz = TRUE,leaflab="none")
#-----------------------Function Library--------------------------------

classify.atoms<-function(x)
{
  A<-0
  B<-0
  loc<-NULL
  for(i in 1:length(x))
  {
    if(x[i] > 50)
    {
      A<-A+1
      loc[i]<-"S"
    }else{
      loc[i]<-"N"
    }
    if(x[i] < 10)
    {
      B<-B+1
    }
  }
  A.pct<-A/length(x)
  B.pct<-B/length(x)
  return(list(A.pct=A.pct,B.pct=B.pct,Surface.loc=loc))
}

test.classify.atoms<-classify.atoms(pdb$atoms$temp)
test.classify.atoms

difference.quotient<- function(x, y) {
  n <- length(x)
  i1 <- 1:2
  i2 <- (n-1):n
  c(diff(y[i1]) / diff(x[i1]), 
    (y[-i1] - y[-i2]) / (x[-i1] - x[-i2]),
    diff(y[i2]) / diff(x[i2]))
}

