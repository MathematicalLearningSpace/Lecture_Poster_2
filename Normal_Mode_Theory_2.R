library(gdata)
library(bio3d)
library(igraph)
library(sna)
library(ips)
library(phangorn)
library(proteomics)
library(dcGOR)
library(MDplot)
library(UniProt.ws)
library(circlize)
library(BioPhysConnectoR)
library(protr)
library(seqinr)
library(Biostrings)
library(Peptides)

library(PearsonDS)
library(xtable)

#------------------------------------------Data------------------------------------------------------

pdb_2VTB<- read.pdb("2VTB.pdb")

#-----------------------------------------Elastic Network Models------------------------------------

pdb_2VTB.open<- trim.pdb(pdb_2VTB, 
                         atom.select(pdb_2VTB, chain="A"))

modes.ANM.2VTB.A<- nma(pdb_2VTB.open, ff="anm", 
                       mass=FALSE, temp=NULL)
modes.calpha.2VTB.A<- nma(pdb_2VTB.open, ff="calpha", mass=FALSE, temp=NULL)
modes.pfanm.2VTB.A<- nma(pdb_2VTB.open, ff="pfanm", mass=FALSE, temp=NULL)
modes.reach.2VTB.A<- nma(pdb_2VTB.open, ff="reach", mass=FALSE, temp=NULL)
modes.sdenm.2VTB.A<- nma(pdb_2VTB.open, ff="sdenm", mass=FALSE, temp=NULL)

select.2VTB <- atom.select(pdb_2VTB, chain='A', elety='CA')
xyz <- pdb_2VTB$xyz[select.2VTB$xyz]
hessian <- build.hessian(xyz, energy.Function)
modes_2VTB.A2 <- eigen(hessian)

modes_2VTB.A$frequencies
#---------------------------------------Trajectory of Modes------------------------------------------ 

Trajectory.2VTB.mode.7<-mktrj(modes.ANM.2VTB.A, mode=7)
Trajectory.2VTB.mode.8<-mktrj(modes.ANM.2VTB.A, mode=8)
Trajectory.2VTB.mode.9<-mktrj(modes.ANM.2VTB.A, mode=9)
Trajectory.2VTB.mode.10<-mktrj(modes.ANM.2VTB.A, mode=10)

correlation.map<- dccm(modes.ANM.2VTB.A)

desmodes_2VTB.A2 <- density(modes_2VTB.A2$vectors)
desmodes_2VTB.A<-density(modes_2VTB.A$frequencies)
desmodes_2VTB.D<-density(modes_2VTB.D$frequencies)
desmodes_2VTB.B<-density(modes_2VTB.B$frequencies)

#-----------------------------------------Deformation Energy----------------------------------------

energy.2VTB <- deformation.nma(modes_2VTB.A) 
fluctuations.2VTB <- fluct.nma(modes_2VTB.A, mode.inds=seq(7,11))
fluctuations.2VTB.reach <- fluct.nma(modes.reach.2VTB.A, mode.inds=seq(7,11))
#-----------------------------------------Tables----------------------------------------------------




#----------------------------------------Figures----------------------------------------------------

Figure.1<-plot(modes.ANM.2VTB.A, sse=pdb_2VTB.open)
Figure.2<-plot(modes.calpha.2VTB.A,sse=pdb_2VTB.open)
Figure.3<-plot(modes.pfanm.2VTB.A,sse=pdb_2VTB.open)
Figure.4<-plot(modes.reach.2VTB.A,sse=pdb_2VTB.open)
Figure.5<-plot(modes.sdenm.2VTB.A,sse=pdb_2VTB.open)

Figure.6<-plot(Trajectory.2VTB.mode.7)
title(main='Mode Trajectories for 2VTB',lty=1,col=1)
lines(Trajectory.2VTB.mode.8,lty=2,col=2)
lines(Trajectory.2VTB.mode.9,lty=3,col=3)
lines(Trajectory.2VTB.mode.10,lty=4,col=4)
cols <- c("black","red", "blue","green")
legend("topright", legend=c("Mode=7",
                            "Mode=8", "Mode=9",
                            "Mode=10")
       , bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:4)

Figure.7<-plot(correlation.map, sse=pdb_2VTB.open, contour=F, 
               col.regions=bwr.colors(20), at=seq(-1,1,0.1) )

Figure.8<-hist(modes_2VTB.A$frequencies)

xlim <- range(0,1)
ylim <- range(0,desmodes_2VTB.A$y)
c2VTBDCol <- rgb(0,1,0,0.2)
c2VTBACol <- rgb(0,0,1,0.2)
c2VTBBCol <- rgb(1,0,0,0.2)
Figure.9<-plot(desmodes_2VTB.A, xlim = xlim, ylim = ylim, xlab = 'Frequencies',
     main = 'Distribution', 
     panel.first = grid())
polygon(desmodes_2VTB.A, density = -1, col = c2VTBACol)
polygon(desmodes_2VTB.D, density = -1, col = c2VTBDCol)
polygon(desmodes_2VTB.B, density = -1, col = c2VTBBCol)
legend('topleft',c('2VTB Chain A','2VTB Chain D','2VTB Chain B'),
       fill = c(c2VTBACol,c2VTBDCol,c2VTBBCol), bty = 'n',
       border = NA)

Figure.9<-plot(fluctuations.2VTB,lty=1,col=1)
lines(fluctuations.2VTB.reach,lty=2,col=2)
title(main='Fluctuations for 2VTB')
cols <- c("black","red")
legend("topright", legend=c("Mode:7-11",
                            "Mode:7-11 Reach")
       , bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:4)
#---------------------------------------Function Library--------------------------------------------
energy.Function<-function(x)
{
    #Theoretical 
    ifelse(x > 10,0,1)

}

