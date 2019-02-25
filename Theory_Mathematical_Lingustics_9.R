#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(RPostgreSQL);library(RSQLite);library(RSQLServer);library(RODBC);require(stringr);require(stringi);library(XML);library(tm);library(plyr)
library(xtable);library(lipsum);library(wordcloud)

#----------------------------------------Data--------------------------------------------------------------------
q <- c(1,2)
#----------------------------------------Search Engine Design------------------------------------------
gene.journal.articles<-lapply(q, function(x)  article.search(terms = x))

#--------------------Student Formatted Notes from Independent Readings of Journal Articles-----------------------
Journal.Articles.Curated <- read_csv("Journal_Articles_Curated.txt", col_names = FALSE)
View(Journal.Articles.Curated)
                              
Journal.Articles.Curated.df<-as.data.frame(Journal.Articles.Curated)
Journal.Articles.Curated.author.list<-JournalArticles.1.df$authors
Journal.Articles.Curated.Title.list<-JournalArticles.1.df$titles
Journal.Articles.Curated.journal.list<-JournalArticles.1.df$journals
Journal.Articles.Curated.references.df<-data.frame()
#--------------------------------------Data Processing-----------------------------------------------------------

gene.corpus <- Corpus(VectorSource(gene.journal.articles))
(dtm <- DocumentTermMatrix(gene.corpus))
journal.articles.content.Frequency<-findFreqTerms(dtm, lowfreq=1, highfreq = 5)

journal.articles.content.Frequency.2<-findFreqTerms.modification(dtm, lowfreq=1, highfreq = 5)
journal.articles.content.Frequency.2

tdm.1 <- TermDocumentMatrix(gene.corpus, control = list(wordLengths = c(3,10)))
gene.corpus.associations<-findAssocs(tdm.1,"vitae",.9)

tdm.1.m <- as.matrix(tdm.1)
tdm.1.m.sort <- sort(rowSums(tdm.1.m),decreasing=TRUE)
tdm.1.m.sort.df <- data.frame(word = names(tdm.1.m.sort),freq=tdm.1.m.sort)
row.names(tdm.1.m.sort.df)<-tdm.1.m.sort.df$word
article.table<-table(tdm.1.m.sort.df$freq)

#---------------------------------------Tables-------------------------------------------------------------------

Table.1<-xtable(Journal.Articles.df)
Table.2<-xtable(as.data.frame(journal.articles.content.Frequency.2$rowNorm))
Table.3<-xtable(tdm.1.m.sort.df)

#---------------------------------------Figures------------------------------------------------------------------

Figure.1<-wordcloud(tdm.1.m.sort.df$word,tdm.1.m.sort.df$freq, scale=c(8,.2),
                    min.freq=2,max.words=Inf, 
                    random.order=FALSE, rot.per=.15)
Figure.2<-hist(journal.articles.content.Frequency.2$rowNorm)

#---------------------------------------Function Library---------------------------------------------------------
article.search<-function(terms)
{
  url = ""
  args<-list(terms = terms)
  #out <- GET(url, query = args, ...)
  res<-lipsum[terms]
  return(res)
}
findFreqTerms.modification<-function (x, lowfreq = 0, highfreq = Inf) 
{
  x <- t(x)
  rs <- slam::row_sums(x)
  rn <-slam::row_norms(x)
  list(rowSum=names(rs[rs >= lowfreq & rs <= highfreq]),rowNorm=rn)
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
                              
                              
