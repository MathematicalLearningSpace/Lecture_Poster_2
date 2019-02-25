#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(europepmc);library(rentrez);library(tm);library(slam);library(NLP);library(openNLP);library(XML);library(RTextTools);library(Rstem)
library(topicmodels);library(wordcloud);library(wordnet);library(proxy);library(plyr);library(wordmatch);library(xml2);library(XML);library(xtable)
library(igraph);library(netgen);library(ggplot2);library(ergm);library(qdap);library(reutils);library(ontologyIndex);library(RNeXML)
library(Rstem);library(RODBC)

#-----------------------------------Data-------------------------------------------------------------------------
abstracts<-c("This is a math paper 1.",
             "This is a math paper 2.")
#----------------------------------Process Data---------------------------------------------------------
document<-abstracts
RDCorpus <- Corpus(VectorSource(document));
RDCorpus = tm_map(RDCorpus, content_transformer(tolower));
RDCorpus = tm_map(RDCorpus, content_transformer(removePunctuation));
RDCorpus = tm_map(RDCorpus, content_transformer(removeNumbers));
RDStopwords = c(stopwords('english'), "available", "via");
idx = which(RDStopwords == "r");
RDStopwords = RDStopwords[-idx];
RDCorpus = tm_map(RDCorpus, content_transformer(removeWords), RDStopwords);
dictCorpus = RDCorpus;
RDCorpus = tm_map(RDCorpus, stemDocument);
RDDtm = DocumentTermMatrix(RDCorpus, control = list(minWordLength = 3));
WordList<-c("math","paper","avian")
WordListTDM <- TermDocumentMatrix(RDCorpus, control = list(dictionary = WordList))
k=5
rownames(RDDtm)
freq <- colSums(as.matrix(RDDtm))
length(freq)
ord <- order(freq,decreasing=TRUE)
freq[ord]
burnin <- 4000
iter <- 2000
thin <- 500
seed <-list(2003,5,63,100001,765)
nstart <- 5
best <- TRUE

#----------------------------------Latent Dirchlet Models---------------------------------------------------------
lda.solution <-LDA(RDDtm,k, method='Gibbs', 
             control=list(nstart=nstart, seed = seed, best=best, 
                          burnin = burnin, iter = iter, thin=thin))
lda.solution.posterior <- posterior(lda.solution,RDDtm)
ldaOut.topics <- as.matrix(topics(lda.solution))

lda.solution.terms <- as.matrix(terms(lda.solution,6))
topicProbabilities <- as.data.frame(lda.solution@gamma)
topic1ToTopic2 <- lapply(1:nrow(RDDtm),function(x)
  sort(topicProbabilities[x,])[k]/sort(topicProbabilities[x,])[k-1])
topic2ToTopic3 <- lapply(1:nrow(RDDtm),function(x)
  sort(topicProbabilities[x,])[k-1]/sort(topicProbabilities[x,])[k-2])
#----------------------------------Correlated Topic Models---------------------------------------------------------
ctm.solution <- CTM(RDDtm, k = 3)
lambda<-0.5
ctm.matrix.1<-build_graph(ctm.solution, lambda, and = TRUE)
ctm.matrix.2<-build_graph(ctm.solution, lambda, and = FALSE)
perplexity(ctm.solution)
ctm.solution@gamma

#-----------------------------------Tables-------------------------------------------------------------------------
Table.1<-xtable(as.data.frame(ldaOut.terms))
Table.2<-xtable(topicProbabilities)
Table.3<-xtable(as.data.frame(topic1ToTopic2))
Table.4<-xtable(as.data.frame(topic2ToTopic3))

#----------------------------------Figures------------------------------------------------------------------------
Figure.1<-hist(lda.solution.posterior$terms)
Figure.2<-hist(lda.solution.posterior$topics)
Figure.3<-hist(topicProbabilities$V1)
Figure.4<-hist(ctm.solution@gamma)

#---------------------------------Function Library--------------------------------------------------------------
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
