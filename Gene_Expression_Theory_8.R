#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-------------------------------------------------------R API ---------------------------------------------
library(HDMD);library(xtable);library(Matrix);library(GSEABase);library(Peptides);library(ActiveDriver)
library(Bio3d);library(protr);library(seqinr);library(CHNOSZ);library(gss);library(PearsonDS);library(fitdistrplus)
library(igraph);library(sets);library(stringr);library(stringi)
#--------------------------------------------Data-phosphorylated proteins and phospholipids ------------------------------

data(Sachs) #protein expression levels in human immune system cells under stimulations.
Sachs.names.description<-c("praf",	"Raf phosphorylated at S259",
"pmek",	"Mek1/mek2 phosphorylated at S217/S221.",
"plcg",	"Phosphorylation of phospholipase C-?? on Y783.",
"pip2",	"Phophatidylinositol 4,5-biphosphate.",
"pip3",	"	Phophatidylinositol 3,4,5-triphosphate.",
"p44.42","Erk1/erk2 phosphorylated at T202/Y204.",
"pakts473","AKT phosphorylated at S473.",
"pka","Phosphorylation of of protein kinase A substrates on 3 sites.",
"pkc","Phosphorylation of of protein kinase C substrates on S660.",
"p38","Erk1/erk2 phosphorylated at T180/Y182.",
"pjnk","Erk1/erk2 phosphorylated at T183/Y185.")

View(Sachs)
AminoAcids<-c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AA.seq.1<-sample(AminoAcids,10)
AA.seq.2<-sample(AminoAcids,10)
AA.seq.3<-sample(AminoAcids,10)

AA.sequence.group.matrix.1<-matrix(c(AA.seq.1,
                              AA.seq.2,
                              AA.seq.3),nrow=3,ncol=10)

#-------------------------------------------Analysis-------------------------------------------
Notes<-c("Flow cytometry measure protein expression level and measures of protein-modification states-phosphorylation")

summary(Sachs)

pearson.MLE.praf <- pearsonFitML(Sachs$praf)
pearson.MLE.pmek <- pearsonFitML(Sachs$pmek)
pearson.MLE.plcg <- pearsonFitML(Sachs$plcg)
pearson.MLE.pip2 <- pearsonFitML(Sachs$pip2)
pearson.MLE.pip3 <- pearsonFitML(Sachs$pip3)
pearson.MLE.p44.42 <- pearsonFitML(Sachs$p44.42)
pearson.MLE.pakts473 <- pearsonFitML(Sachs$pakts473)
pearson.MLE.pka <- pearsonFitML(Sachs$pka)
pearson.MLE.pkc <- pearsonFitML(Sachs$pkc)
pearson.MLE.p38 <- pearsonFitML(Sachs$p38)
pearson.MLE.pjnk <- pearsonFitML(Sachs$pjnk)

pearson.N<-length(Sachs$p38)
pIpars <- list(a=pearson.MLE.pkc$a, b=pearson.MLE.pkc$b, location=pearson.MLE.pkc$location, scale=pearson.MLE.pkc$scale) 
pIVpars <- list(m=pearson.MLE.p38$m, nu=pearson.MLE.p38$nu, location=pearson.MLE.p38$location, scale=pearson.MLE.p38$scale)
error.pearson.1<-rpearsonI(pearson.N,params=pIpars)
error.pearson.4<-rpearsonIV(pearson.N,params=pIVpars)

Pearson.data.df<-data.frame()
Pearson.data.df<-rbind(c(pearson.MLE.praf$type),
                       c(pearson.MLE.pmek$type),
                       c(pearson.MLE.plcg$type),
                       c(pearson.MLE.pip2$type),
                       c(pearson.MLE.pip3$type),
                       c(pearson.MLE.p44.42$type),
                       c(pearson.MLE.pakts473$type),
                       c(pearson.MLE.pka$type),
                       c(pearson.MLE.pkc$type),
                       c(pearson.MLE.p38$type),
                       c(pearson.MLE.pjnk$type))
colnames(Pearson.data.df)<-c("Distribution")
rownames(Pearson.data.df)<-c(Sachs.names.description[seq(1,22,2)])

distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c(Sachs.names.description[1],pearson.MLE.praf$type,empMoments.Version.1(Sachs$praf)),
                                 c(Sachs.names.description[3],pearson.MLE.pmek$type,empMoments.Version.1(Sachs$pmek)),
                                 c(Sachs.names.description[5],pearson.MLE.plcg$type,empMoments.Version.1(Sachs$plcg)),
                                 c(Sachs.names.description[7],pearson.MLE.pip2$type,empMoments.Version.1(Sachs$pip2)),
                                 c(Sachs.names.description[9],pearson.MLE.pip3$type,empMoments.Version.1(Sachs$pip3)),
                                 c(Sachs.names.description[11],pearson.MLE.p44.42$type,empMoments.Version.1(Sachs$p44.42)),
                                 c(Sachs.names.description[13],pearson.MLE.pakts473$type,empMoments.Version.1(Sachs$pakts473)),
                                 c(Sachs.names.description[15],pearson.MLE.pka$type,empMoments.Version.1(Sachs$pka)),
                                 c(Sachs.names.description[17],pearson.MLE.pkc$type,empMoments.Version.1(Sachs$pkc)),
                                 c(Sachs.names.description[19],pearson.MLE.p38$type,empMoments.Version.1(Sachs$p38)),
                                 c(Sachs.names.description[21],pearson.MLE.pjnk$type,empMoments.Version.1(Sachs$pjnk))
                                 
)
colnames(distribution.moments.1.df)<-c("Variable","Type", "Moment 1", "Moment 2","Moment 3","Moment 4","Theta")

AA.seq.1.Molecular.Entropy<-MolecularEntropy(AA.seq.1,"AA")
AA.seq.2.Molecular.Entropy<-MolecularEntropy(AA.seq.2,"AA")
AA.seq.3.Molecular.Entropy<-MolecularEntropy(AA.seq.3,"AA")

AA.seq.1.MolecularMI<-MolecularMI(AA.seq.1,"AA")
AA.sequence.group.matrix.1.Molecular.MI<-MolecularMI(AA.sequence.group.matrix.1,"AA")

AA.sequence.Properties.df<-data.frame()
AA.sequence.Properties.df<-rbind(c(AA.seq.1.Molecular.Entropy$H),
                                 c(AA.seq.2.Molecular.Entropy$H),
                                 c(AA.seq.3.Molecular.Entropy$H))
colnames(AA.sequence.Properties.df)<-c("Entropy")

#---------------------------------------------Graph--------------------------------------------

graph.V<-c(Sachs.names.description[seq(1,22,2)])
graph.1<-make_full_graph(length(graph.V))

E(graph.1)$weight <- seq_len(ecount(graph.1))

graph.1.edges<-E(graph.1)
graph.2<-delete_edges(graph.1,seq(1, length(graph.1.edges), by = 2))
graph.3<-delete_edges(graph.1,seq(1, length(graph.1.edges), by = 1))
graph.4<-delete_edges(graph.1,seq(1, length(graph.1.edges), by = 3))

E(graph.2)$weight<-seq(0,1,length.out=ecount(graph.2))
E(graph.3)$weight<-seq(0,1,length.out=ecount(graph.3))
E(graph.4)$weight<-rep(1,ecount(graph.4))

graph.0<-make_empty_graph(length(graph.V))

graph.0.1<-add_edges(graph.0,c(1,1))
graph.0.2<-add_edges(graph.0,c(1,2))
graph.0.3<-add_edges(graph.0,c(1,3))

E(graph.0.1)$weight<-seq(0,1,length.out=ecount(graph.0.1))
E(graph.0.2)$weight<-seq(0,1,length.out=ecount(graph.0.2))
E(graph.0.3)$weight<-seq(0,1,length.out=ecount(graph.0.3))

#---------------------------------------------Tables-------------------------------------------
Table.1<-xtable(AA.sequence.Properties.df)
Table.2<-xtable(Pearson.data.df)
Table.3<-xtable(distribution.moments.1.df)

#----------Figures for Presentation in the Classroom------------------------------------------

Figure.1<-plot(diag(AA.sequence.group.matrix.1.Molecular.MI), 
               type = "h", 
               ylab="Functional Entropy", xlab="site")

par(mfrow = c(3,3), mar = c(0.5,1,0.5,0.5))
Figure.2A<-hist(Sachs$praf)
Figure.2B<-hist(Sachs$pmek)
Figure.2C<-hist(Sachs$plcg)
Figure.2D<-hist(Sachs$pip2)
Figure.2E<-hist(Sachs$pip3)
Figure.2F<-hist(Sachs$p44.42)
Figure.2G<-hist(Sachs$pakts473)
Figure.2H<-hist(Sachs$pka)
Figure.2I<-hist(Sachs$pkc)
Figure.2J<-hist(Sachs$p38)
Figure.2K<-hist(Sachs$pjnk)

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.3<-plot(graph.0)
Figure.4<-plot(graph.0.1)
Figure.5<-plot(graph.0.2)
Figure.6<-plot(graph.0.3)

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.7<-plot(graph.1)
Figure.8<-plot(graph.2)
Figure.9<-plot(graph.3)
Figure.10<-plot(graph.4)

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.11<-plot_dendrogram(cluster_fast_greedy(graph.1))
Figure.12<-plot_dendrogram(cluster_fast_greedy(graph.2))
Figure.13<-plot_dendrogram(cluster_fast_greedy(graph.3))
Figure.14<-plot_dendrogram(cluster_fast_greedy(graph.4))

Figure.15<-plot_dendrogram(cluster_walktrap(graph.1,steps=5))

#---------------------------------------------References--------------------------------------

Reference.1<-c("Sachs, K., Perez, O., Pe'er, D., Lauffenburger, D. A., and Nolan, G. P. (2005)", 
"Causal protein-signaling networks derived from multiparameter single-cell data.", 
"Science, 308 (5732), 523-529.")

#----------------Function Library to be modified by Students-------------------------------
empMoments.table.caption<-function(nouns, references)
{
  verbs<-c("")
  objects<-c("")
  sentences<-NULL
  for(i in 1:length(nouns))
  {
    sentences[i]<-stri_join(nouns[i]," ",verbs[1]," ",objects[1], " from ",references[i]," .")
  }
  return(sentences)
}

empMoments.Version.1<-function (x) 
{
  n <- length(x)
  moment.1 <- mean(x)
  moment.2  <- var(x) * (n - 1)/n
  if (moment.2 > 0) 
    moment.3 <- sum((x - moment.1)^3/sqrt(moment.2)^3)/n
  else moment.3 <- 0
  if (moment.2 > 0) 
    moment.4 <- sum((x - moment.1)^4/moment.2^2)/n
  else moment.4 <- 0
  if(moment.4> 0)
    theta.1<-moment.3/moment.4 
  else theta.1<-0
  c(average = moment.1, 
    variance = moment.2, 
    skewness = moment.3, 
    kurtosis = moment.4, 
    theta=theta.1)
}

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
