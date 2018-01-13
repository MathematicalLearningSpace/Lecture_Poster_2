library(WikipediR)
library(xtable)
library(RWeka)
library(tm)
library(RCurl)
library(topicmodels)
library(wordcloud)
library(stringi)
library(corpora)
library(zipfR)
#-----------------------------------------Data-----------------------------------------------------------

optimization.data<-getURLContent('https://cran.r-project.org/web/views/Optimization.html')
optimization<-remove.HTML.markup(optimization.data)
article.translate <- as.String(optimization)
article.Sentences<-strsplit(article.translate, "[.]")
article.Words = lapply(article.translate, strsplit, "[[:space:]]")

article.content.1 <-page_content("en","wikipedia", page_name = "Mathematical optimization")
article.content.2 <-page_content("en","wikipedia", page_name = "genetic algorithm")

article.1<-sapply(article.content.1,remove.HTML.markup)
article.2<-sapply(article.content.2,remove.HTML.markup)

biomolecules.list<-remove.HTML.markup(getURLContent('https://en.wikipedia.org/wiki/List_of_biomolecules'))
chem.term<-grep("c+",biomolecules.list, perl=TRUE, value=TRUE)

Nucleolus<-getURLContent('https://en.wikipedia.org/wiki/Nucleolus')

html.files <- list.files(pattern="\\.(htm|html)$")

eukaryotic.Cell<-c(Nucleolus,Nucleus,Ribosome,Vesicle,Roughendoplasmicreticulum,
                   Golgiapparatus,Cytoskeleton, Smoothendoplasmicreticulum, Mitochondrion,
                   Vacuole, Cytosol,Lysosome,Centrosome,CellMembrane)
eukaryotic.Cell.df<-as.data.frame(eukaryotic.Cell)


save(eukaryotic.Cell.df,file='eukaryotic.Cell.df.RData')

#-----------------------------------------Data Processing------------------------------------------------

docs<-c(article.1,article.2)
docs<-VectorSource(docs)
corpus <- VCorpus(docs)
corpus<- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("english"))
corpus <- tm_map(corpus, stripWhitespace)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, PlainTextDocument)

dtm <- DocumentTermMatrix(corpus)
freqs_dtms <- colSums(as.matrix(dtm))
freqs_dtms
freqs_dtms[head(order(freqs_dtms,decreasing = TRUE),75)]
freqs_dtms[head(order(freqs_dtms,decreasing = FALSE),75)]

findFreqTerms(dtm, lowfreq=500)
findFreqTerms(dtm, highfreq=500)

tdm.1 <- TermDocumentMatrix(corpus, control = list(wordLengths = c(3,10)))
findAssocs(tdm.1,"optimization",.9)

articles.cluster<-kmeans(tdm.1, 4)

tdm.1.m <- as.matrix(tdm.1)
tdm.1.m.sort <- sort(rowSums(tdm.1.m),decreasing=TRUE)
tdm.1.m.sort.df <- data.frame(word = names(tdm.1.m.sort),freq=tdm.1.m.sort)
row.names(tdm.1.m.sort.df)<-tdm.1.m.sort.df$word
article.table<-table(tdm.1.m.sort.df$freq)

#-----------------------------------------Tables--------------------------------------------------------

Table.1<-xtable(tdm.1.m.sort.df)

#----------------------------------------Figures-------------------------------------------------------

Figure.1<-wordcloud(tdm.1.m.sort.df$word,tdm.1.m.sort.df$freq, scale=c(8,.2),
                    min.freq=3,max.words=Inf, 
                    random.order=FALSE, rot.per=.15)


#----------------------------------------Function Library----------------------------------------------

remove.HTML.markup<-function(s)
  { 
    doc.1<-NULL
    doc.1 <- htmlTreeParse(paste("<!DOCTYPE html>", s),asText = TRUE, trim = FALSE)
    xmlValue(xmlRoot(doc.1))
  } 
