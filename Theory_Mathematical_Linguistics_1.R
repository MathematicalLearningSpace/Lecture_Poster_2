#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------------------R API ----------------------------------------------------------
library(tm);library(stringi);library(stringr);library(utils);library(CHNOSZ);library(PearsonDS);library(xtable)
library(bio3d);library(Peptides);library(jsonlite);library(rjson);library(xtable);library(seqinr)
#--------------------------------------------Data------------------------------------------------------------------------
data(thermo)
data(aa.table)
#---------------------------------Example Test Sequence for Classroom Discussion--------------------------------
test.seq<-c('MAANLSRNGPALQEAYVRVVTEKSPTDWALFTYEGNSNDIRVAGTGEGGLEEMVEELNSGK
            VMYAFCRVKDPNSGLPKFVLINWTGEGVNDVRKGACASHVSTMASFLKGAHVTINARAEEDVEPECI
            MEKVAKASGANYSFHKESGRFQDVGPQAPVGSVYQKTNAVSEIKRVGKDSFWAKAEKEEENRRLEE
            KRRAEEAQRQLEQERRERELREAARREQRYQEQGGEASPQSRTWEQQQEVVSRNRNEQGSTCASL
            QESAVHPREIFKQKERAMSTTSISSPQPGKLRSPFLQKQLTQPETHFGREPAAAISRPRADL
            PAEEPAPSTPPCLVQAEEEAVYEEPPEQETFYEQPPLVQQQGAGSEHIDHHIQGQGLSGQGL
            CARALYDYQAADDTEISFDPENLITGIEVIDEGWWRGYGPDGHFGMFPANYVELIE')
test.seq.2<-str_to_upper(s2c(test.seq))
test.seq.pattern<-stri_locate_all_regex("ACAGAGACTTTAGATAGAGAAGA", "(?=AGA)")

article.project.files <- list.files(patt='*.*pdf$')
#----------------------------Examples of WDRs---------------------------------------------
keywords<-c('coronin','WD40','WDR5','WDR64')
#-----------------Formatted Student Notes and Examples for the Classroom-------
WDR5.df<- read.delim("WDR5_HUMAN.txt")
WDR5<-read.fasta('WDR5_HUMAN.fasta')
motif = c("VM....CI")
motif.find(motif, test.seq)
#-----------------------------------------Language Processing for Student Work Here------------------------------------------------------------
#-----------------------------------------WD40 Example---------------------------------------------------------------------------
#-----------------------------------------Beta strands D - WD40 beta proteins-------------------------------------------
WD40.tetrad<- 'DHSW' 
WD40.Strand.D.1.protein.sequence<-'ALKFTL'
WD40.Strand.D.2.protein.sequence<-'KFEKTI'
WD40.Strand.D.3.protein.sequence<-'KCLKTL'
WD40.Strand.D.4.protein.sequence<-'KCLKTL'
WD40.Strand.D.5.protein.sequence<-'QCLKTL'
WD40.Strand.D.6.protein.sequence<-'KCLKTY'
WD40.Strand.D.7.protein.sequence<-'EIVQKL'

aa.test <- seq2aa('test',WD40.Strand.D.1.protein.sequence)
pf.test<- protein.formula(aa.test)
cf.test<-as.chemical.formula(pf.test)

pf.test.df<-data.frame(pf.test)

Hydrogen.Carbon.ratio <- function(d) d$H/d$C
Oxygen.Carbon.ratio <- function(d) d$O/d$C
Nitrogen.Carbon.ratio <- function(d) d$N/d$C
Sulfur.Carbon.ratio <- function(d) d$S/d$C

Hydrogen.Carbon.ratio(pf.test.df)

#-----------------------------------------WDR5---------------------------------------------------------------------------

sentence.1<- stri_paste("There are ", length(WDR5.df$Repeats)," repeats.")
aa.strand.d.Face.Side <- seq2aa('WDR5.Strand.D',WDR5.df$Strand_d)
aa.strand.c.Core.Inner.Interior.1 <- seq2aa('WDR5.Strand.c',WDR5.df$Strand_c)
aa.strand.b.Core.Inner.Interior.2<- seq2aa('WDR5.Strand.b',WDR5.df$Strand_b)
aa.strand.a.Core.Inner <- seq2aa('WDR5.Strand.a',WDR5.df$Strand_a)

#----------------------------------------Average composition per beta strand---------------------------------------------
aa.composition.average.strand.df<-data.frame()
aa.composition.average.strand.df<-rbind(aasum(aa.strand.d.Face.Side, average=TRUE),
                                        aasum(aa.strand.c.Core.Inner.Interior.1, average=TRUE),
                                        aasum(aa.strand.b.Core.Inner.Interior.2, average=TRUE),
                                        aasum(aa.strand.a.Core.Inner, average=TRUE))
aa.loop.da.Face.Top <- seq2aa('WDR5.loop.da',WDR5.df$Loop_da)
aa.loop.bc <- seq2aa('WDR5.loop.bc',WDR5.df$Loop_bc)
aa.loop.ab.Face.Bottom.1 <- seq2aa('WDR5.loop.ab',WDR5.df$Loop_ab)
aa.loop.cd.Face.Bottom.2 <- seq2aa('WDR5.loop.cd',WDR5.df$Loop_cd)
#---------------------------------------Average composition per beta strand----------------------------------------------
aa.composition.average.loop.df<-data.frame()
aa.composition.average.loop.df<-rbind(aasum(aa.loop.da.Face.Top, average=TRUE),
                                      aasum(aa.loop.bc, average=TRUE),
                                      aasum(aa.loop.ab.Face.Bottom.1, average=TRUE),
                                      aasum(aa.loop.cd.Face.Bottom.2, average=TRUE))

WDR.components<-c('Bulge','Face.Top','Face.Side','Face.Bottom','Core')
WDR.components.df<-data.frame()
WDR.components.df<-rbind(aa.composition.average.strand.df,
                         aa.composition.average.loop.df)

#---------------------------------------Tables-------------------------------------------------------------------------
Table.1<-xtable(aa.composition.average.loop.df)
Table.2<-xtable(WDR.components.df)
#---------------------------------------Figures-----------------------------------------------------------------------
Figure.1<-aa_Stats<-AAstat(test.seq.2)
title("Residue Statistics for Test Sequence")
Figure.2<-barplot(aa_Stats$Compo, xlab="Amino Acid",
                  ylab = "Frequency",
                  col="lightcyan",
                  angle=85,
                  legend.text = c("Test Sequence"),
                  args.legend = list(x = "topright"))
title(main=paste("AA Distribution for \n","Test Sequence"))

#---------------------------------------Function Library--------------------------------------------------------------
carbon.Oxidation.average<-function(Hn,Nn,On,Sn,Cn,molecule.charge=0)
{
  coa<-(-Hn+3*Nn+2*On+2*S*n+molecule.charge)/Cn
  return(coa)
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
