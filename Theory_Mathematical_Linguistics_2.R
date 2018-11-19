#------------------------------------------------R API --------------------------------------------
library(europepmc);library(rentrez);library(tm);library(slam);library(NLP);library(openNLP);library(XML)
library(RTextTools);library(Rstem);library(topicmodels);library(wordcloud);library(wordnet);library(proxy)
library(plyr);library(wordmatch);library(xml2);library(XML);library(xtable);library(igraph);library(netgen)
library(ggplot2);library(ergm);library(qdap);library(reutils);library(ontologyIndex);library(RNeXML);library(Rstem)
library(RODBC)

#---------------------Data Sets for Students in the Classroom----------------------------------------------------------------
StatementOfResearchHypothesis<-c("Avian Spatial Intelligence")

channel<-odbcConnect("AvianDB",rows_at_time = 1)
SQLSelect<-"SELECT [AvianName]
FROM [AvianDB].[dbo].[AvianMathPapers]
Where [AvianName] Like '%ZebraFinch%'"
ZebraFinch.DF<- as.data.frame(sqlQuery(channel,paste(SQLSelect)))
odbcClose(channel)
View(ZebraFinch.DF)

#------------------------------------------------Update Database-----------------------------------------------------
abstracts<-c("This is an avian math paper 1.",
             "This is an avian math paper 2.",
             "This is an avian math paper 3.",
             "This is an avian paper on spatial intelligence.",
             "This is an avian paper on the cerebellum.")
abstractWords = lapply(abstracts, strsplit, "[[:space:]]") 
abstractWords = lapply(abstractWords, function(x) wordStem[[1]])
lapply(abstractWords, function(x) x[x %in% stopWords])
length(abstracts)
AbstractsArray<-as.array(abstracts)
abstracts.DF=as.data.frame(AbstractsArray)
TempAbstractPreDataMathIQ<-abstracts.DF
View(abstracts.DF)

channel<-odbcConnect("ScientificPublicationsDB",rows_at_time = 1)
sqlSave(channel, TempAbstractPreDataMathIQ, rownames = FALSE)
sqlQuery(channel, 'INSERT INTO [dbo].[Abstract]
         ([AbstractContextPre])
         select AbstractsArray from SciencePublicationsDB.dbo.TempAbstractPreData')
odbcClose(channel)
abstracts.DF<-as.data.frame(sqlQuery(channel, 'Select AbstractContextPre from SciencePublicationsDB.dbo.Abstract'))

#--------------------------------------------Process the Abstracts----------------------------------------------------
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
RDCorpus = tm_map(RDCorpus,stemDocument);
#RDCorpus = tm_map(RDCorpus, stemCompletion, dictionary=dictCorpus);

AbstractPostArray<-as.data.frame(as.matrix(RDCorpus))
View(AbstractPostArray)

channel<-odbcConnect("ScientificPublicationsDB",rows_at_time = 1)
sqlSave(channel, TempAbstractPreData, rownames = FALSE)
sqlQuery(channel, 'INSERT INTO [dbo].[Abstract]
         ([AbstractContextPost])
         select AbstractPostArray from ScientificPublicationsDB.dbo.TempAbstractPostData')
odbcClose(channel)

#-------------------------------------------Term Matrix---------------------------------------------------------

RDDtm = DocumentTermMatrix(RDCorpus, control = list(minWordLength = 3));
RDTdm <- TermDocumentMatrix(RDCorpus)
WordList<-c("math","paper","avian")
WordListTDM <- TermDocumentMatrix(RDCorpus, control = list(dictionary = WordList))
termDocMatrix=inspect(RDTdm)
termDocMatrix[termDocMatrix>=1] <- 1

termMatrix <- termDocMatrix %*% t(termDocMatrix)
desired <- as.data.frame(as.table(termMatrix))
termMatrix<-termMatrix[1:4,1:4]
#-------------------------------------------Rank the Matrix------------------------------------------------------
m <- as.matrix(RDTdm)
v <- sort(rowSums(m), decreasing=TRUE)
N<-5
RankedTermDocMatrix<-head(v, N)
RankedTermDocMatrix[RankedTermDocMatrix>=1] <- 1
RankedTermMatrix <- RankedTermDocMatrix %*% t(RankedTermDocMatrix)

research.keywords<-c("avian","math","paper")
RDAssociations<-findAssocs(WordListTDM,research.keywords, c(0.9,0.9,0.9))
RDAssociationDF<-data.frame()
#RDAssociationDF<-head(as.data.frame(findAssocs(RDTdm,research.keywords,0.1)),n=5)
RDAssociationDF<-cbind(RDAssociationDF,rownames(RDAssociationDF))
View(RDAssociationDF)

#------------------------------------------------Graph Theory--------------------------------------------------------
g <- graph.data.frame(RDAssociationDF, directed = TRUE)
g <- graph.adjacency(termMatrix, weighted=T, mode = 'undirected')
gRanked <- graph.adjacency(RankedTermMatrix, weighted=T, mode = 'undirected')
#-------------------------------------------------Tables--------------------------------------------------------------
Table.1<-xtable(RDAssociationDF)
#-------------------------------------------------Figures to be added in the classroom-------------------------------------------------------------
Figure.1<-plot(g,layout=layout.kamada.kawai, vertex.size=1, vertex.color="green")

#------------------------------------------------Function Library----------------------------------------------------

