#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(wordcloud);library(VarfromPDB);library(tm);library(mRMRe);library(europepmc);library(rentrez);library(aRxiv)
#--------------------------Data--------------------------------------
data(cgps)
cgps.ge
cgps.ic50
#------------------------------Student Formatted Notes from Independent Readings---------------------------
publications <- read_csv("search-results/publications.csv")
TERT.publications.df<-as.data.frame(publications[grep("TERT", publications$Abstract),])
TERT.publications.df$Title
TERT.publications.df$Abstract

chromosome.study<-TERT.publications.df[grep("chromosome",TERT.publications.df$Abstract),]
chromosome.5<-cgps.annot[grep("^5",cgps.annot$chromosome_name),]

chromosome.5$band

Gene.study.df<-data.frame()
Gene.study.df<-cbind(chromosome.5$hgnc_symbol,chromosome.5$description)
colnames(Gene.study.df)<-c("Gene Name","Description")

gene.expression.df<-as.data.frame(cgps.ge)
x<-list()
#-----------------------------------------------SQL for First Search-----------------
for(i in 1:length(chromosome.5))
{
  x[i]<-gene.expression.df[rownames(chromosome.5)[i]]
}
y<-list()
for(i in 1:length(chromosome.5$hgnc_symbol))
{
  y[i]<-entrez_search(db="pubmed",term=paste(chromosome.5$hgnc_symbol[i]," AND (2017[PDAT])"),retmax = 10)
}
#--------------------------Parse XML Data Structure---------------------------
search.pubmed <- entrez_fetch(db = "pubmed", id = y[14][[1]],rettype = "xml", parsed = T)

abstracts = xpathApply(search.pubmed, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))
names(abstracts) <- y[14][[1]]
literature.review.references<-y[14][[1]]
Literature.Review<-filter_abstract(literature.review.references,abstracts)

#----------------------Read and Translate Articles and Abstracts---------------
corp <- Corpus(VectorSource(abstracts))
corp <- tm_map(corp, removePunctuation)
corp <- tm_map(corp, content_transformer(tolower))
corp <- tm_map(corp, removeNumbers)
corp <- tm_map(corp, function(x)removeWords(x,stopwords()))
corp <- tm_map(corp, function(x)removeWords(x,c("author",
                                                "manuscript",
                                                "pubmed","pmc",
                                                "page","available")))
term.matrix <- TermDocumentMatrix(corp)
findAssocs(term.matrix,chromosome.5$hgnc_symbol[14],.5)
term.matrix <- as.matrix(term.matrix)
Top.Terms<-head(sort(rowSums(term.matrix),decreasing = TRUE))
Top.Terms

#-------Tables For Descriptive Statistics and Testing Differences in Moments-------

Table.1<-xtable(as.data.frame(moment.analytics(x,chromosome.5$hgnc_symbol,FALSE)))

#-----------------------Figures-----------------------------------------------------
Figure.1<-plot(chromosome.5$jetset.overall, xaxt = "n",type="l") axis(1, 1:length(chromosome.5$hgnc_symbol),chromosome.5$hgnc_symbol)

#-------------------------------Figure Group 1-----------------------------------
par(mfcol = c(3, 3))
for(i in 1:9)
{
  hist(x[i][[1]], xlab="Expression",
       main = paste("Histogram of",chromosome.5$hgnc_symbol[i] ))
  mtext(outer = FALSE, side = 3, paste("mean=",mean(x[i][[1]]),
                                       "sd=",sd(x[i][[1]])), cex = 0.5)
}
#-------------------------------Figure Group 2-----------------------------------
par(mfcol = c(3, 3))
for(i in 10:18)
{
  hist(x[i][[1]], xlab="Expression",
       main = paste("Histogram of",chromosome.5$hgnc_symbol[i] ))
  mtext(outer = FALSE, side = 3, paste("mean=",mean(x[i][[1]]),
                                       "sd=",sd(x[i][[1]])), cex = 0.5)
}
#-------------------------------Figure Group 3-----------------------------------
par(mfcol = c(2, 2))
Figure.2<-wordcloud(corp, min.freq=2,random.order=FALSE)
Figure.3<-wordcloud(TERT.publications.df$Abstract, min.freq=2,random.order=FALSE)
Figure.4<-comparison.cloud(term.matrix,max.words=2500,random.order=FALSE)
Figure.5<-commonality.cloud(term.matrix,max.words=2500,random.order=FALSE)
#----------------------References----------------------------------

Reference.1<-c("Stephen Bridgett, James Campbell, Christopher J. Lord, and Colm J. Ryan,",
"CancerGD: A Resource for Identifying and Interpreting Genetic Dependencies in Cancer",
"Cell Syst. 2017 Jul 26; 5(1): 82-86.e3. doi:  10.1016/j.cels.2017.06.002")

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

#----------------------Function Library----------------------------
moment.analytics<-function(x,x.names,plot.x=FALSE)
{
  moment.1<-NULL
  for(i in 1:length(x))
  {
    moment.1[i]<-mean(x[i][[1]])
    
  }
  if(plot.x==TRUE)
    par(mfcol = c(1, 1))
    plot(moment.1,type="l",xaxt = "n",xlab="Gene",ylab="Expression")
    axis(1, 1:length(x.names),x.names,cex=0.25)
  
  return(moment.1)
}

test.moment<-moment.analytics(x,chromosome.5$hgnc_symbol,TRUE)
test.moment

filter_abstract<-function(x,y)
{
  LiteratureReview<-NULL
  sNbr<-length(x)
  for(i in 1:length(x))
  {
    s<- as.character(y[i])
    s<-strsplit(s,split = "[\\.?!] ")
    s<-as.data.frame(s)
    nbrSentences<-nrow(s)
    if(nbrSentences > 6)
    {
      LiteratureReview[i]<-c(paste(x[i]," states ",
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
