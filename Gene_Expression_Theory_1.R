#-----------------------------------R API ----------------------------------------------------
library(fractal); library(fracdiff);library(timeSeries);library(RandomFields);library(boot);library(vars)
library(urca);library(PearsonDS);library(stringr);library(stringi);library(corrplot);library(ggplot2);library(xtable)

pearson.N<-100
#-------------------------------------Moment Parameters for the Distributions-------------------------------------------------

p0pars <- list(mean=1, sd=1)
pIpars <- list(a=1, b=1, location=1, scale=1) 
pIIpars <- list(a=1, location=1, scale=1) 
pIIIpars <- list(shape=1, location=1, scale=1)
pIVpars <- list(m=3, nu=1, location=1, scale=1)
pVpars <- list(shape=1, location=1, scale=1) 
pVIpars <- list(a=1, b=1, location=1, scale=1) 
pVIIpars <- list(df=10, location=1, scale=1)

#-------------------------------------Generate Random variables from the Distributions----------------------------------------

error.pearson.0<-rpearson0(pearson.N,params=p0pars)
error.pearson.1<-rpearsonI(pearson.N,params=pIpars)
error.pearson.2<-rpearsonII(pearson.N,params=pIIpars)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
error.pearson.4<-rpearsonIV(pearson.N,params=pIVpars)
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
error.pearson.6<-rpearsonVI(pearson.N,params=pVIpars)
error.pearson.7<-rpearsonVII(pearson.N,params=pVIIpars)

#-------------------------------Create Gene Theoretical Data Set based on CDF Hypothesis 
#--------------for Analysis Workflow based on Simulated Network Data for the Classroom-------------------------------------

MOV10<-window(error.pearson.0,start = 25,end=75)
CSNK2B<-window(error.pearson.1,start = 25,end=75)
AMMECR1<-window(error.pearson.2,start = 25,end=75)
VRK1<-window(error.pearson.3,start = 25,end=75)
BRCA1<-window(error.pearson.4,start = 25,end=75)
P53<-window(error.pearson.5,start = 25,end=75)
GNE<-window(error.pearson.6,start = 25,end=75)
GABPA<-window(error.pearson.7,start = 25,end=75)

Gene.Study<-as.matrix(cbind(MOV10,CSNK2B,AMMECR1,VRK1,
                            BRCA1,P53,GNE,GABPA))
Gene.Study.Names<-c("MOV10","CSNK2B","AMMECR1","VRK1","BRCA1","P53","GNE","GABPA")

#-------------------------------Five different expression profiles (G1/S, S, G2, G2/M, and M/G1)-----

CellCycle.Type<-c('G1/S','S','G2','G2/M','M/G1')
CellCycle.Type.G1<-c('E2F1')
CellCycle.Type.S<-c('RFC4')
CellCycle.Type.G2<-c('CDC2')
CellCycle.Type.M<-c('STK15')
CellCycle.Type.MG1<-c('PTTG1')

Gene.Expression.Table.df<-data.frame()
Gene.Expression.Table.df<-cbind(CellCycle.Type.G1,CellCycle.Type.S,
                                CellCycle.Type.G2,CellCycle.Type.M,
                                CellCycle.Type.MG1)
colnames(Gene.Expression.Table.df)<-c(CellCycle.Type[1],
                                      CellCycle.Type[2],
                                      CellCycle.Type[3],
                                      CellCycle.Type[4],
                                      CellCycle.Type[5])

Gene.Matrix.Correlation<-cor(Gene.Study)
Gene.Matrix.Correlation.Test<-correlation.Test(Gene.Study,0.99)

#--------------------------------------Stationarity Testing----------------------------------------
x<-Gene.Study[,1]
lc.df <- ur.df(y=x, lags=3, type='trend')
ers.Test <- ur.ers(x, type="DF-GLS", model="const", lag.max=4)
pp.Test <- ur.pp(x, type="Z-tau", model="trend", lags="short")
sp.Test <- ur.sp(x, type="tau", pol.deg=1, signif=0.01)
za.Test <- ur.za(x, model="both", lag=2)

NameOfTest<-c(summary(ers.Test)@test.name,summary(pp.Test)@test.name,
              summary(lc.df)@test.name, summary(sp.Test)@test.name,
              summary(za.Test)@test.name)

TestStatisic<-c(summary(ers.Test)@teststat,summary(pp.Test)@teststat,
                summary(lc.df)@teststat, summary(sp.Test)@teststat,
                summary(za.Test)@teststat)  

CriticalValue1PCT<-c(summary(ers.Test)@cval[1], summary(pp.Test)@cval[1],
                     summary(lc.df)@cval[1],summary(sp.Test)@cval[1],
                     summary(za.Test)@cval[1])

Test.Stationarity<-NULL
Test.Stationarity<-rbind(Test.Stationarity,NameOfTest)
Test.Stationarity<-rbind(Test.Stationarity,TestStatisic)
Test.Stationarity<-rbind(Test.Stationarity,CriticalValue1PCT)

#--------------------------------------Hurst Estimation for Gene Expression-------------------------

Hurst.methods <- c("aggabs","aggvar","diffvar","higuchi")
Hurst.methods.spectral <- c("standard","smoothed","robinson")
Gene.Hurst<-c('Gene','Higuchi',"Standard","Smoothed","Robinson")
Hurst.Gene<-NULL
Hurst.Gene<-rbind(Hurst.Gene,Gene.Hurst)
Hurst.Statistic.Gene<-NULL
for (i in 1:length(Gene.Study.Names)){
  walk<-cumsum(Gene.Study[,i])
  z<- lapply(Hurst.methods[4], function(method, walk){
    hurstBlock(ifelse1(method=="higuchi",diff(walk),walk), method=method)
  },walk=walk )
  names(z) <- Hurst.methods[4]
  Hurst.Spec.z <- lapply(Hurst.methods.spectral, function(method, walk){
    hurstSpec(walk, method=method, sdf.method="multitaper")
  },walk=walk )
  names(Hurst.Spec.z) <- Hurst.methods.spectral
  Hurst.Statistic.Gene<-c(Gene.Study.Names[i],format(z$higuchi[1], digits=4),
                          format(Hurst.Spec.z$standard[1], digits=4),
                          format(Hurst.Spec.z$smoothed[1], digits=4),
                          format(Hurst.Spec.z$robinson[1], digits=4))
  Hurst.Gene<-rbind(Hurst.Gene,Hurst.Statistic.Gene)
}
#------------------------------Tables-------------------------------------------------------------------------

Table.1<-xtable(Gene.Expression.Table.df)
Table.2<-xtable(Test.Stationarity)
Table.3<-xtable(Hurst.Gene)
#-------------Figures for Presentation in the Classroom-----------------------------------------------------------------------
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", "cyan", "#007FFF", "blue","#00007F"))

Figure.1<-plot(Gene.Study[,1],type='l', xlab="Hours",ylab='Expression')
title(main='Gene Expression Cell Cycle',lty=1,col=1)
lines(Gene.Study[,2],lty=2,col=2)
lines(Gene.Study[,3], lty=3,col=3)
lines(Gene.Study[,4], lty=4,col=4)
lines(Gene.Study[,5], lty=5,col=5)
lines(Gene.Study[,6], lty=6,col=6)
lines(Gene.Study[,7], lty=7,col=7)
cols <- c("black","red","blue","green","yellow","orange","gray","brown","purple")
legend("topright", legend=c(Gene.Study.Names[1],
                            Gene.Study.Names[2],
                            Gene.Study.Names[3],
                            Gene.Study.Names[4],
                            Gene.Study.Names[5],
                            Gene.Study.Names[6],
                            Gene.Study.Names[7]), bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:7)

Figure.2<-corrplot(Gene.Matrix.Correlation, 
                   order="hclust", addrect=3,col=col1(100),t1.cex=0.8)
title(main="Gene Correlation Plot")
#-----------------------------Function Library for Modification in the Classroom-----------------------------------
correlation.Test <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      test <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- test$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- test$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- test$conf.int[2]
    }
  }
  return(list(p.mat, 
              lowCI.mat, 
              uppCI.mat,
              conf.level)
         )
}
