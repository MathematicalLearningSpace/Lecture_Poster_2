library(gsubfn)
library(zipfR)
library(koRpus)
library(KoNLP)
library(textreuse)
library(boilerpipeR)
library(RCurl)
library(XML)
library(tm)
library(rebus)
library(stringi)
library(arules)
library(qdapRegex)
library(qdap)
library(hash)
library(corpora)
library(textcat)
library(textir)
library(RTextTools)
library(kernlab)
library(lda)
library(topicmodels)
library(lsa)
library(readability)
library(syllable)
library(sylcount)

#-----------------------------------------------------Data------------------------------------------------------------------
abstract.test<-c("This is an abstract.","This is the definition of a circulant graph.")
Definitions.Data <- data.frame(
  text = c("Definition 1. This is a definition of a circulant graph.", 
           "Definition 2. This is a definition of a graph.", 
           "Definition 3. This is a defintion of eigenvalues."),
  group = c("A", "A", "B")
)
math.dictionary<-design.dict(abstract.test[1])

poster.features <- list(
  sentences=18, 
  words=556, 
  letters=c(all=2918, l1=19, l2=92, l3=74, l4=80, l5=51, l6=49), 
  syllables=c(all=974, s1=316, s2=116), 
  punct=78, 
  all.chars=3553, 
  prepositions=74, 
  conjunctions=18, 
  pronouns=9, 
  foreign=0, 
  TTR=0.5269784, 
  Bormuth.NOL=192, 
  Dale.Chall.NOL=192, 
  Harris.Jacobson.NOL=240, 
  Spache.NOL=240)

#----------------------------------------------Word Associations-----------------------------------------------------
corpus.doc<-Corpus(VectorSource(abstract.test))
tdm <- TermDocumentMatrix(corpus.doc)
words <- c("graph", "circulant")
corr <- c(0.7, 0.75)
AssocTerms <- findAssocs(tdm, words, corr)
ListAssocTerms <- lapply(AssocTerms, function(x) data.frame(terms = names(x),
                                                            cor = x, stringsAsFactors=FALSE))
ListAssocTerms[1:10]$graph

poster.readability<-readability.num(poster.features, index="all")

readability(c(abstract.test[1], abstract.test[2]))
readability(paste0(Definitions.Data$text[1], Definitions.Data$text[1], collapse=" "))

readability_word_stats(abstract.test)
readability_word_stats_by(Definitions.Data$text, Definitions.Data$group)
#----------------------------------------------------Tables----------------------------------------------------------------

#----------------------------------------------------Figures---------------------------------------------------------------

#----------------------------------------------------Function Library------------------------------------------------------
design.dict <- function(text, max.vocab=10000) {
  text <- strsplit(text, '')
  diction1 <- list()
  indx <- 1
  for (c in text[[1]]) {
    if (!(c %in% names(diction1))) {
      diction1[[c]] <- indx
      indx <- indx + 1
    }
  }
  return (diction1)
}