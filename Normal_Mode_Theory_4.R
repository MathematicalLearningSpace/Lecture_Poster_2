library(bio3d)
library(xtable)
library(Peptides)
library(stringi)
library(xtable)
#------------------------------------------Data---------------------------------------------------------

WDR<-c(WDR1, WDR5,WDR10, WDR12, WDR13, WDR16, 
       WDR17,WDR18, WDR19, WDR20, WDR21A, 
       WDR21C, WDR22, WDR23, WDR24, WDR25, 
       WDR26, WDR27, WDR3, WDR31, WDR32, WDR33, 
       WDR34, WDR35, WDR36, WDR37, WDR38, WDR4, 
       WDR40A, WDR40B, WDR40C, WDR41, WDR42A, WDR42B, WDR43, 
       WDR44, WDR46, WDR47, WDR48, WDR49, WDR5, WDR51A, WDR51B, 
       WDR52, WDR53, WDR54, WDR55, WDR57, WDR59, WDR5B, 
       WDR6, WDR60, WDR61, WDR62, WDR63, WDR64, 
       WDR65, WDR66, WDR67, WDR68, WDR69, WDR7, 
       WDR70, WDR72, WDR73, WDR74, WDR75, WDR76, 
       WDR77, WDR78, WDR79, WDR8, WDR81, WDR82, 
       WDR85, WDR86, WDR88, WDR89, WDR90, WDR91,WDR92)

WDR5.pdb<-read.pdb("WDR5_HUMAN_1.pdb")
COR1A.pdb<-read.pdb("COR1A_Human_1.pdb")
NormalMode.Study<-c(WDR5.pdb,COR1A.pdb)

#----------------------------------------Load the force fields-----------------------------------
ff.1 <- load.enmff('calpha')
ff.2 <- load.enmff('anm')
ff.3 <- load.enmff('pfanm')
ff.4 <- load.enmff('sdenm')
ff.5 <- load.enmff('reach')

ca.inds <- atom.select(WDR5.pdb, 'calpha')
xyz <- WDR5.pdb$xyz[ca.inds$xyz]
dists <- dm.xyz(xyz, mask.lower=FALSE)
force.matrix <- apply(dists, 1, ff.1.mutation)

#----------------------------------------Build Hessian Matrix-------------------------------------
WDR5.hessian<-build.hessian(WDR5.pdb$xyz[ ca.inds$xyz ],pfc.fun=ff.1)
WDR5.hessian.1 <- build.hessian(WDR5.pdb$xyz[ ca.inds$xyz ], 
                     pfc.fun=ff.1.mutation)

#----------------------------------------Spectral Properties------------------------------------------
modes.hessian.base <- eigen(WDR5.hessian, symmetric=TRUE)
modes.hessian.1 <- eigen(WDR5.hessian.1, symmetric=TRUE)

#-----------------------------------------Mode Robustness---------------------------------------------
a <- nma(WDR5.pdb, ff="calpha") 
b <- nma(WDR5.pdb, ff="anm")
c <- nma(COR1A.pdb, ff="calpha") 
d <- nma(COR1A.pdb, ff="anm")
#-----------------------------------------Root mean square inner product (RMSIP)-----------------------
modes.RMSIP.WDR5 <- rmsip(a, b)
modes.RMSIP.COR1A <- rmsip(c, d)

#------------------------------------------Tables-------------------------------------------------------


#-----------------------------------------Figures-------------------------------------------------------
r<-seq(0,10,by=0.1)
par(mfcol = c(2, 3))
Figure.1<-plot(r,ff.1(r), col="Blue")
Figure.2<-plot(r,ff.1.mutation(r),col="red")
Figure.3<-plot(r,ff.2(r),col="black")
Figure.4<-plot(r,ff.3(r),col="green")
Figure.4A<-plot(r,ff.1.mutation.2(r,0.5),col="purple")

Figure.5<-plot(modes.hessian.base$values,,lty=1,col=1)
title(main='Eigenvalues for WDR5 Force Field Mutation')
lines(modes.hessian.1$values,lty=2,col=2)
cols <- c("black","red")
legend("topright", legend=c("FF Base Eigenvalues",
                            "FF Mutation Eigenvalues")
       , bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:2)

#----------------------------------------Function Library-----------------------------------------------

ff.1.mutation<-function (r, rmin = 2.9, ...) 
{
  if (!is.null(rmin)) 
    r[(r < rmin)] = rmin
  a <- 128 * 10^4
  b <- 8.6 * 10^2
  c <- 2.39 * 10
  ifelse(r < 4, b * (r) - c, a * (r)^(-6))
}
ff.1.mutation.2<-function(r,p)
{
  ifelse( r>8, 0, r^(-p) )
}