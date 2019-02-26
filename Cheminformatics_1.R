#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------------R API -----------------------------------------------------

library(CHNOSZ);library(rcdk);library(leaps);library(caret);library(igraph);library(GA);library(webchem)
library(markovchain);library(mlr);library(mlrMBO);library(mlrCPO);library(rpubchem);library(klaR)
library(data.table);library(fuzzyjoin);library(tidytext);library(dplyr);library(tm)
library(pROC);library(MASS);library(rpubchem);library(readr);library(rjson);library(e1071);library(ROCR)







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
