#---------------------------------------------------R API --------------------------------------------------------
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);library(PearsonDS);library(xtable)

#------------------------------Data-------------------------------------------------------

pdb_2VTB<- read.pdb("2VTB.pdb")
pdb_2VTB.open<- trim.pdb(pdb_2VTB, 
                         atom.select(pdb_2VTB, chain="A"))
modes.calpha.2VTB.A<- nma(pdb_2VTB.open, ff="calpha", mass=FALSE, temp=NULL)
correlation.map<- dccm(modes.calpha.2VTB.A)

#----------------------------------Sequence------------------------------------------------
aa_2VTB.seq<-pdbseq(pdb_2VTB)
motif = c("G....GKS")
motif.find(motif, aa_2VTB.seq)
biounit <- biounit(pdb_2VTB)
biounit

#-----------------------------------Network-----------------------------------------------

network.1<-cna(correlation.map, cutoff.cij=0.9)
network.2<-cna(correlation.map, cluster.method='walk', cutoff.cij=0.9)
network.3<-cna(correlation.map, cluster.method ='greed', cutoff.cij=0.9)

network.table.1<-summary(network.1)
attributes(network.1)
network.communties.table.1<-table( network.1$communities$members )
network.communties.table.2<-table( network.2$communities$members )
network.communties.table.3<-table( network.3$communities$members )

network.1.max<-max(network.1$communities$modularity)
network.1.tree <- community.tree(network.1, rescale=TRUE) 
network.1.pruned <- prune.cna(network.1, size.min = 30)

#-----------------------------------Path in the Network---------------------

network.1.path <- cnapath(network.1, from=14, to=25, k=10)

#-----------------------------------Graph Theory---------------------------

adjM_Full<-as_adjacency_matrix(network.1$network)
adjM_Communities<-as_adjacency_matrix(network.1$community.network)
adjM_Pruned<-as_adjacency_matrix(network.1.pruned$community.network)

g1 <- graph_from_adjacency_matrix(adjM_Full)
g2 <- graph_from_adjacency_matrix(adjM_Communities)
g3 <- graph_from_adjacency_matrix(adjM_Pruned)

#-----------------------------Tables-----------------------------------------------------
Table.1<-xtable(network.table.1$tbl)

#----------------------------Figures to be Presented in the Classroom-----------------------------------------------------

par(mfcol = c(2, 2), mar = c(0, 0, 0, 0)) 
Figure.1<-plot(network.1, pdb_2VTB.open, full = TRUE,vertex.size=5, 
     vertex.label.cex=0.7) 
Figure.2<-plot(network.2, pdb_2VTB.open)
Figure.3<-plot(network.3, pdb_2VTB.open)
Figure.4<-plot(network.1.tree$num.of.comms,
      network.1.tree$modularity, 
      xlab="Communities", 
      ylab="Modularity")
Figure.5<-plot(network.1.pruned,pdb_2VTB.open)
Figure.6<-print(network.1.path, plot=TRUE)
lo<-layout_in_circle(g1)
Figure.7<-plot(g1,layout=layout_with_kk, vertex.color="Blue",vertex.size=3,
     vertex.label.dist=1,edge.arrow.size=0.2,main='', 
     vertex.label.cex=0.6,margin=c(-0.05,-0.05,-0.05,0))
Figure.8<-plot(g2,layout=lo, vertex.color="Red",vertex.size=5,
     vertex.label.dist=1,edge.arrow.size=0.2,main='', 
     vertex.label.cex=0.6,margin=c(-0.15,-0.15,-0.15,0))
Figure.9<-plot(g3,layout=lo, vertex.color="Green",vertex.size=5,
     vertex.label.dist=1,edge.arrow.size=0.2,main='', 
     vertex.label.cex=0.6,margin=c(-0.15,-0.15,-0.15,0))
Figure.10<-hist(network.communties.table.3)
#----------------------------Function Library-------------------------------------------


