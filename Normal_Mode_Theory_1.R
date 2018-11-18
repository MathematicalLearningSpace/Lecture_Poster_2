#----------------------------------------------------R API------------------------------------
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);library(PearsonDS);library(xtable)

#-----Example Data from Protein Database Retrieved by Students in the Classroom------------------------------

pdb_2VTB<- read.pdb("2VTB.pdb")

centerMass_2VTB<-com(pdb_2VTB)
aa_2VTB.seq<-pdbseq(pdb_2VTB)

tor <- torsion.pdb(pdb_2VTB)

inds <- atom.select(pdb_2VTB, "calpha")
a.inds_2VTB <- atom.select(pdb_2VTB, chain="A")

#------------------Modification of the Hydrophobic Moment by Students--------------------------------------------
sequence.angle.100<-hmoment(seq = c2s(pdb_2VTB$seqres), angle = 100, window = 11)
sequence.angle.160<-hmoment(seq = c2s(), angle = 160, window = 11)

AA.Index.Table.df<-data.frame()
AA.Index.Table.df<-cbind(sequence.angle.100,sequence.angle.160)
colnames(AA.Index.Table.df)<-c("HMoment 100","HMoment 160")
#-------------------------------------------------Contact Map---------------------------------------------
ref.cont <- cmap(pdb_2VTB$xyz[inds$xyz], dcut=6, scut=3 )
g_Contact<-graph_from_adjacency_matrix(ref.cont)
#-------------------------------------------------Factors-------------------------------------------------
bfactor.2VTB<-pdb_2VTB$atom$b[pdb_2VTB$calpha]
bfactor<-data.frame(cbind(bfactor.2VTB))
#-------------------------------------------------Density------------------------------------------------
dens.2VTB <- density(bfactor[,1])
#------------------------------------------------Normal Mode Analysis-------------------------------------
pdb_2VTB.A.open<- trim.pdb(pdb_2VTB, 
                         atom.select(pdb_2VTB, chain="A"))
pdb_2VTB.B.open<- trim.pdb(pdb_2VTB, 
                           atom.select(pdb_2VTB, chain="B"))
pdb_2VTB.D.open<- trim.pdb(pdb_2VTB, 
                           atom.select(pdb_2VTB, chain="D"))

modes_2VTB.A <- nma(pdb_2VTB.A.open)
modes_2VTB.B <- nma(pdb_2VTB.B.open)
modes_2VTB.D <- nma(pdb_2VTB.D.open)

nv_modes.7_2VTB.A<- normalize.vector(modes_2VTB.A$modes[,7])
nv_modes.7_2VTB.B<- normalize.vector(modes_2VTB.B$modes[,7])
nv_modes.7_2VTB.D<- normalize.vector(modes_2VTB.D$modes[,7])
#------------------------------------------------Tables----------------------------------------------------------

Table.1<-xtable(AA.Index.Table.df)

#--------------Figures for Presentation and Discussion--------------------------------------------------------
Figure.1<-plot(g_Contact,layout=layout_with_kk, vertex.color="Blue",vertex.size=3,
     vertex.label.dist=1,edge.arrow.size=0.2,main='', 
     vertex.label.cex=0.3,margin=c(-0.35,-0.35,-0.35,0))

Figure.2<-plot.cmap(ref.cont)

Figure.3<-plot(tor$phi, tor$psi, xlab='phi',ylab='psi',
               main='Ramachandran Plot 2VTB',col='Blue')
legend("topright", legend=c("2VTB")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)

Figure.4<-hist(bfactor[,1], xlim=c(0,80), col= rgb(1,0,0,0.2))


c2VTBCol <- rgb(1,0,0,0.2)
Figure.5<-plot(dens.2VTB, xlab = 'B-Factor',
     main = 'Distribution of B-Factor for 2VTB Paralogs', 
     panel.first = grid())
polygon(dens.2VTB, density = -1, col = c2VTBCol)
legend('topleft',c('2VTB'),
       fill = c(c2VTBCol), bty = 'n',
       border = NA)

Figure.6<-plot.bio3d(pdb_2VTB$atom$b[pdb_2VTB$calpha], sse=pdb_2VTB, 
           typ="l",ylim=c(25,80),ylab="B-factor")
title(main=' Residue Temperature Factors for 2VTB')

Figure.7<-plot(nv_modes.7_2VTB.A, ylab="Normalized Mode 7",xlab='Residue')
points(nv_modes.7_2VTB.D,col=rgb(1,0,0,0.8))
points(nv_modes.7_2VTB.B,col=rgb(0,0,1,0.8))
## add a legend in the corner
legend('bottomleft',c('2VTB Chain A','2VTB Chain D','2VTB Chain B'),
       fill = c('black','red','blue'), bty = 'n',
       border = NA)
#-----------------------------------------------Function Library------------------------------------------------
