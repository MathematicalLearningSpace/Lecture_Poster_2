#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(githubinstall);library(devtools);library(europepmc);library(rentrez);library(aRxiv);library(XML);library(pubmed.mineR);library(stringi)
library(stringr);library(xtable);library(NLP);library(openNLP);library(tm);library(wordcloud)

#----------------------------------------------------Data--------------------------------------------------------------
P53.data<- epmc_search(query = 'P53',limit = 250)

View(P53.data)
#----------------------------------------------------100 most cited atricles from PLOS ONE-----------------------------
PLOS.data<-epmc_search(query = 'ISSN:	1932-6203', sort = 'CITED desc')

#----------------------------------------------------arXiv Data for Mathematical Biology--------------------------------------------------------
Equations.Differential<- arxiv_search("Differential Equations", sort_by="updated",ascending=FALSE)
DE.Abstracts<-Equations.Differential$abstract
math.biology.count<-arxiv_count("92Bxx")
arXiv.categories<-arxiv_cats

#----------------------------------------------------R Methods for Math Biology----------------------------------------
Pubmed.R.Methods<- entrez_search(db="pubmed", term="R Language", retmax=100)
Pubmed.R.Methods.Equations.Differential.Stochastic<- entrez_search(db="pubmed", term="Stochastic Differential Equations", retmax=100)
Pubmed.R.Methods.Equations.Differential.Delayed<- entrez_search(db="pubmed", term="Delayed Differential Equations", retmax=100)
Pubmed.R.Methods.Equations.Differential.Fractional<- entrez_search(db="pubmed", term="Fractional Differential Equations", retmax=100)

#----------------------------------------------------Pubmed------------------------------------------------------------
year <- 2000:2017
papers <- sapply(year, search_PubMed_year, term="Cell Cycle G1 S Signaling", USE.NAMES=FALSE)
paper<-epmc_ftxt(P53.data$pmcid[2])
paper <- xmlParse(paper)
poster.references <- c( P53.data$pmid[1],P53.data$pmid[2],P53.data$pmid[3])

search.pubmed <- entrez_fetch(db = "pubmed", id = poster.references,
                             rettype = "xml", parsed = T)

#----------------------------------------------------Citations---------------------------------------------------------

Authors<-xpathApply(search.pubmed, '//PubmedArticle//Article//AuthorList', function(x)
  xmlValue(xmlChildren(x)$Author))
Author<-xpathApply(search.pubmed, '//PubmedArticle//Article//AuthorList/Author', function(x)
  xmlValue(xmlChildren(x)$LastName))
#Author<-lapply(Author, strsplit, "[[:space:]]") 
TitleOfPaper<-xpathApply(search.pubmed, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$ArticleTitle))
JournalTitle<-xpathApply(search.pubmed, '//PubmedArticle//Article//Journal', function(x)
  xmlValue(xmlChildren(x)$Title))
JournalIssue<-xpathApply(search.pubmed, '//PubmedArticle//Article//Journal', function(x)
  xmlValue(xmlChildren(x)$JournalIssue))
#Citation for the Introduction
citation.Intro<-as.array(Author)
#Citation for the Bibliography
citation.BiB<-paste(Authors,TitleOfPaper,JournalTitle,JournalIssue)

abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))
names(abstracts) <- poster.references

abstracts.df<-data.frame()
abstracts.df<-rbind(abstracts[[1]],
                    abstracts[[2]])
colnames(abstracts.df)<-c("P53")

#-----------------------------------------------Generate Draft of Literature Review------------------------------

LiteratureReview_P53<-filter_abstract(poster.references,abstracts)
wordList<-generate_Word_list(poster.references,abstracts)
trigrams<-ngrams(wordList, 3L)

#---------------------------------------------Analyze Abstracts------------------------------
corp <- Corpus(VectorSource(abstracts))
corp <- tm_map(corp, removePunctuation)
corp <- tm_map(corp, content_transformer(tolower))
corp <- tm_map(corp, removeNumbers)
corp <- tm_map(corp, function(x)removeWords(x,stopwords()))

term.matrix <- TermDocumentMatrix(corp)
term.matrix <- as.matrix(term.matrix)
colnames(term.matrix) <- poster.references

write.csv(LiteratureReview_P53, file = "LiteratureReview_P53.csv")

#----------------------------------------------------Tables-------------------------------------------------------------

Table.1<-xtable(as.data.frame(LiteratureReview_P53))

#-----------------------------------------------------Figures-------------------------------------------------------------
par(mfcol = c(2, 2))
Figure.1 <-plot(year, papers, type='b', main="Cell Cycle G1 S Signaling")
Figure.2<-wordcloud(abstracts$`29323737`,min.freq=1, random.color=TRUE,ordered.colors=TRUE)
Figure.3<-wordcloud(abstracts$`29244835`,min.freq=1,random.color=TRUE,ordered.colors=TRUE)
Figure.4<-wordcloud(abstracts$`29311662`,min.freq=1,random.color=TRUE,ordered.colors=TRUE)
par(mfcol = c(1, 1))
Figure.5<-comparison.cloud(term.matrix,max.words=40,random.order=FALSE)
Figure.6<-commonality.cloud(term.matrix,max.words=40,random.order=FALSE)

#-------------------------------------------------Function Library----------------------------------------------------

search_Pubmed_year <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="pubmed", term=query, retmax=0)$count
}
generate_Word_list<-function(x,y)
{
  wordList<-NULL
  sNbr<-length(x)
  for(i in 1:length(x))
  {
    s<- as.character(y[i])
    #Separate the sentences from the paragraph
    w <- strsplit(s, " ", fixed = TRUE)[[1L]]
    ngrams(w, 3L)
  }
  return(w)
}

filter_abstract<-function(x,y)
{
  LiteratureReview<-NULL
  sNbr<-length(x)
  for(i in 1:length(x))
  {
    s<- as.character(y[i])
    #Separate the sentences from the paragraph
    s<-strsplit(s,split = "[\\.?!] ")
    s<-as.data.frame(s)
    nbrSentences<-nrow(s)
    if(nbrSentences > 6)
    {
      LiteratureReview[i]<-c(paste(poster.references[i]," states ",
                                   s[nbrSentences-5,1],
                                   s[nbrSentences-4,1],
                                   s[nbrSentences-3,1]))
    }
  }
  return(LiteratureReview)
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
