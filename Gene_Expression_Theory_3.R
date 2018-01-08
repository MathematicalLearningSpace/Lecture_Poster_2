library(circlize)
library(PearsonDS)
library(igraph)
library(Matrix)
library(matrixStats)
library(matrixcalc)
library(MatrixModels)
library(mathgraph)
library(stringi)
library(deSolve)
library(xtable)
library(Deriv)
library(utils)
library(smoothHR)
library(pspline)
library(bigsplines)
library(formatR)
library(splines)

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
#----------------------------------------------Theoretical Data Set---------------------------------------------------

set.seed(999) 
A = data.frame(factor = sample(letters[1:8], pearson.N, replace = TRUE), 
               x = error.pearson.1, 
               y = error.pearson.3)

data.distribution<-generateRandomBed(nr = pearson.N, 
                                            fun = function(k) error.pearson.6<-rpearsonVI(k,params=pVIpars))

#--------------------------------------------Spline Models---------------------------------------------------------

Spline.linear <- bigspline(A$x,A$y,type="lin")
Spline.cubic <- bigspline(A$x,A$y,type="cub0")

Spline.Table.df<-data.frame()
Spline.Table.df<-cbind(spline.cubic$info,spline.linear$info)
colnames(Spline.Table.df)<-c("Cubic","Linear")

Spline.Prediction.df<-data.frame(x=seq(0,1,length.out=pearson.N))
Spline.Prediction.df<-predict(Spline.cubic,Spline.Prediction.df,se.fit=TRUE)

#---------------------------------------------Linear Models----------------------------------------------------------

model.lm<-lm(A$y ~A$x)
X<-cbind(1,A$x)
model.matrix.lm.solution<-solve(t(X)%*%X)%*%t(X)%*%A$y
summary(model.lm)
anova(model.lm)
model.lm.prediction<-predict(model.lm)

#---------------------------------------------Restriction on Intervals-----------------------------------------------
A.x.Interval<-(A$x-1.5)
A.x.Interval[A.x.Interval<0]<-1
X<-cbind(1,A$x,A.x.Interval)
model.matrix.lm.solution.restriction<-solve(t(X)%*%X)%*%t(X)%*%A$y
model.lm.restriction<-lm(A$y ~A$x+A.x.Interval)
model.lm.restriction.prediction<-predict(model.lm.restriction)
summary(model.lm.restriction)
anova(model.lm.restriction)

f.polynomial(3)
design.x.derivative.1<-Deriv(f.polynomial(3))
design.x.derivative.2<-Deriv(design.x.derivative.2)

#---------------------------------------------Tables-----------------------------------------------------------------

Table.1<-xtable(Spline.Table.df)


#--------------------------------------------Figures-----------------------------------------------------------------
Figure.1<-plot(A$y,type='l', xlab="Hours",ylab='Expression')
title(main='Cubic Spline Prediction',lty=1,col=1)
lines(Spline.Prediction.df$fit,lty=2,col=2)
lines(Spline.Prediction.df$se.fit,lty=3,col=3)
lines(model.lm.prediction,lty=4,col=4)
lines(model.lm.restriction.prediction,lty=5,col=5)
cols <- c("black","red", "blue")
legend("topright", legend=c("Theoretical Data",
                            "Cubic Spline Prediction", "Standard Error",
                            "LM Model","LM Model With Restriction")
                        , bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:5)

par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram()
Figure.1<-circos.genomicTrackPlotRegion(data.distribution, panel.fun = function(region, value, ...) {
  x = (region[[2]] + region[[1]]) / 2
  y = value[[1]]
  loess.fit = loess(y ~ x)
  loess.predict = predict(loess.fit, x, se = TRUE)
  d1 = c(x, rev(x))
  d2 = c(loess.predict$fit + loess.predict$se.fit, rev(loess.predict$fit - loess.predict$se.fit))
  circos.polygon(d1, d2, col = "#CCCCCC", border = NA)
  circos.points(x, y, pch = 16, cex = 0.5)
  circos.lines(x, loess.predict$fit)
  
}, track.height = 0.1)

par(mar = c(1, 1, 1,1), lwd = 0.1, cex = 0.7) 
circos.par("track.height" = 0.1) 
circos.initialize(factors = A$factor, x = A$x)
Figure.2<-circos.trackPlotRegion(factors = A$factor, y = A$y, 
                                 panel.fun = function(x, y) { circos.axis() }) 
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1) 
circos.text(1, 0.5, "right", sector.index = "a")
col = rep(c("#FF0000", "#00FF00"), 4) 
circos.trackPoints(A$factor, A$x, A$y, col = col, pch = 16, cex = 0.8) 
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4) 
circos.trackHist(A$factor, A$x, bg.col = bgcol, col = NA)
circos.trackPlotRegion(factors = A$factor, y = A$y) 
circos.trackLines(A$factor[1:100], A$x[1:100], A$y[1:100], type = "h")
circos.link("a", 0, "b", 0, h = 0.8)
circos.link("b", 0, "c", 0, h = 0.8)
circos.link("c", 0, "d", 0, h = 0.8)
circos.link("d", 0, "e", 0, h = 0.8)
circos.link("e", 0, "f", 0, h = 0.8)
circos.link("f", 0, "g", 0, h = 0.8)
circos.link("g", 0, "h", 0, h = 0.8)
circos.link("h", 0, "a", 0, h = 0.8)

par(mar = c(1, 1, 1, 1)) 
circos.par("canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "clock.wise" = FALSE, "gap.degree" = 1) 
factors = letters[1:4] 
circos.initialize(factors = factors, xlim = c(0, 1)) 
Figure.3<-circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA) 
circos.updatePlotRegion(sector.index = "a", bg.border = "blue") 
x1 = error.pearson.5
y1 = error.pearson.7 
circos.points(x1, y1, pch = 16, cex = 0.5,col="red") 
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = "purple") 
circos.updatePlotRegion(sector.index = "a", bg.border = "green") 
circos.lines(1:100/100, y1, pch = 16, cex = 0.5) 

#--------------------------------------------Function Library-------------------------------------------------------

f.polynomial<-function(polynomial.max) {
  strPoly<-""
  for(i in 1:polynomial.max)
  {
    if(i<polynomial.max)
      strPoly<-stri_join(strPoly,'x','^',as.character(i),"+")
    else
      strPoly<-stri_join(strPoly,'x','^',as.character(i))
  }
  return(strPoly)
}
