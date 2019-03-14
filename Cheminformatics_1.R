#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------------R API -----------------------------------------------------

library(CHNOSZ);library(rcdk);library(leaps);library(caret);library(igraph);library(GA);library(webchem)
library(markovchain);library(mlr);library(mlrMBO);library(mlrCPO);library(rpubchem);library(klaR)
library(data.table);library(fuzzyjoin);library(tidytext);library(dplyr);library(tm)
library(pROC);library(MASS);library(rpubchem);library(readr);library(rjson);library(e1071);library(ROCR)

#--------------------------------\section{Patent Applications}-------------------------------------------
#----------Formatted Student Notes for Classroom Lectures------------------------------------------------
DAHP.Synthase.patent.df<-as.data.frame(read.csv("DAHP.Synthase.patent.csv"))
three.dehydroquinate.synthase.patent.df<-as.data.frame(read.csv("three.dehydroquinate.synthase.patent.csv"))
three.dehydroquinate.dehydratase.patents.df<-as.data.frame(read.csv("three.dehydroquinate.dehydratase.patents.csv"))
chorismate.synthase.patents.df<-as.data.frame(read.csv("chorismate.synthase.patents.csv"))
EPSP.syhthase.patents.df<-as.data.frame(read.csv("EPSP.synthase.patents.csv"))

#---------Formatted Student Notes for Classroom Lectures-----------------------------------------------
#Fingerprint Tanimoto-based 2-dimensional similarity search. 90 %
ASCORBIC.ACID.literature.df<-as.data.frame(read.csv("ASCORBIC_ACID.csv"))
Vitamin.A.literature.df<-as.data.frame(read.csv("Vitamin_A.csv"))
#Fingerprint Tanimoto-based 2-dimensional similarity search 95 %
Vitamin.D.literature.df<-as.data.frame(read.csv("Vitamin_D.csv"))

section.ids <- "Names and Identifiers"
section.ids <- "Computed Descriptors"
section.props <- "Chemical and Physical Properties"
section.computed <- "Computed Properties"
#--------------------------------------Parallel Programming-------------------------------

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
