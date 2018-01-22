library(Biobase)
library(dnet)
library(xtable)
library(rjson)
library(jsonlite)
library(RJSONIO)
library(rentrez)
library(XML)
library(tm)
library(ape)
library(topicmodels)
library(RODBC)
library(stringr)
library(stringi)
#----------------------------------Data-----------------------------------------------------

Ontology.AT.GOBP<-dRDataLoader(genome='At', ontology='GOBP')
Ontology.AT.GOBP.Name<-Ontology.AT.GOBP$set_info$name
Ontology.AT.GOBP.Photosynthesis<-Ontology.AT.GOBP$set_info$name[878:885]
Ontology.AT.GOBP.DNA.Damage<-Ontology.AT.GOBP$set_info$name[19]
Ontology.AT.GOBP.DNA.ID<-Ontology.AT.GOBP$set_info$setID[19]

Ontology.AT.GOBP.DNA.Damage.Genes<-Ontology.AT.GOBP$gs$`GO:0000077`
Ontology.AT.GOBP.DNA.Damage.Genes.JSON<-toJSON(Ontology.AT.GOBP.DNA.Damage.Genes)

#----------------------------------Query the Genes and Publications-------------------------------------------
entrez_dbs()

linked_seq_ids <- entrez_link(dbfrom="gene", id=Ontology.AT.GOBP.DNA.Damage.Genes, db="nuccore")
linked_transripts <- linked_seq_ids$links$gene_nuccore_refseqrna
all_recs <- entrez_fetch(db="nuccore", id=linked_transripts, rettype="fasta")
nchar(all_recs)

all_the_links <- entrez_link(dbfrom='gene', id=Ontology.AT.GOBP.DNA.Damage.Genes, db='all')
all_the_links$links$gene_protein

Gene.Interactions<-interactions.from.gene(Ontology.AT.GOBP.DNA.Damage.Genes)

paper_links <- entrez_link(dbfrom="pubmed", id=Ontology.AT.GOBP.DNA.Damage.Genes, cmd="llinks")
paper_links$linkouts

fetch.pubmed <- entrez_fetch(db = "pubmed", id = all_the_links$links$gene_pmc,rettype = "xml", parsed = T)
abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))

r_search <- entrez_search(db="pubmed", term="R Language and 2015:2017[PDAT]",retmax=40)
DNA.Damage<-entrez_search(db="pubmed", term="DNA Damage and 2015:2017[PDAT]",retmax=40)

fetch.pmc <- entrez_fetch(db = "pubmed", id = r_search$ids,rettype = "xml", parsed = T)
Title = xpathApply(fetch.pmc, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$ArticleTitle))
Abstract = xpathApply(fetch.pmc, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))

fetch.pmc <- entrez_fetch(db = "pubmed", id = DNA.Damage$ids,rettype = "xml", parsed = T)

#----------------Reference, Abstract and Citation Work for Introduction-------------------------------------
Abstract = xpathApply(fetch.pmc, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))
AbstractsArray<-as.array(Abstract)
citation.Intro<-as.array(Author)
#--------------------------------Find term frequency---------------------------------------------------------
vs <- VectorSource(unlist(Abstract))
orig_docs <- Corpus(vs)
docs <- orig_docs
docs <- tm_map(docs, content_transformer(tolower))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, removeWords, stopwords("english"))
RDStopwords = c("the", "and");
docs <- tm_map(docs, content_transformer(removeWords),RDStopwords)

docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs, PlainTextDocument)

#docs.Filter<-tm_filter(docs,FUN = function(x) any(grep("damage", content(x))))
dtm <- DocumentTermMatrix(docs)
d1 = dim(dtm)
dtms <- removeSparseTerms(dtm, .98)
d1=dim(dtms)
#-----------------------------Adjacency Matrix for term Frequency-------------------------------------
termDocMatrix=inspect(dtms)
termDocMatrix[termDocMatrix>=1] <- 1
termMatrix <- termDocMatrix %*% t(termDocMatrix)
termMatrix<-termMatrix[1:10,1:10]
m <- as.matrix(dtms)
v <- sort(rowSums(m), decreasing=TRUE)
N<-10
RankedTermDocMatrix<-head(v, N)
RankedTermDocMatrix[RankedTermDocMatrix>=1] <- 1
RankedTermMatrix <- RankedTermDocMatrix %*% t(RankedTermDocMatrix)

freqs_dtms <- colSums(as.matrix(dtms))
freqs_dtms[head(order(freqs_dtms,decreasing = TRUE),50)]
print("Words more than 500 times")
findFreqTerms(dtm, lowfreq=500)
options(max.print=200)
print("Words less than 2 times")
findFreqTerms(dtm, highfreq=2)
#----------------Associations based on Research Keywords----------------------------------------------
DNA.Assocations<-findAssocs(dtms, 'dna', 0.4)
m = as.matrix(dtms);
v = sort(colSums(m), decreasing=TRUE);
names.1 = names(v);
k = which(names(v)=="dna");
names.1[k] = "dna";
DNA.df = data.frame(word=names.1, freq=v);
research.keywords<-c("cancer")
RDAssociations<-findAssocs(dtms, research.keywords, c(0.75,0.75))
RDAssociations
RDAssociation.df<-head(as.data.frame(findAssocs(dtms,research.keywords,0.5)),n=10)
rownames(RDAssociation.df)[1]
RDAssociation.df<-cbind(RDAssociation.df,rownames(RDAssociation.df))
View(RDAssociation.df)
#-----------------------------Build the Graphs------------------------------------------
g <- graph.data.frame(RDAssociation.df, directed = TRUE)
g.2 <- graph.adjacency(termMatrix, weighted=T, mode = 'undirected')
gRanked <- graph.adjacency(RankedTermMatrix, weighted=T, mode = 'undirected')
V(g.2)$label <- V(g.2)$name
V(gRanked)$label <- V(gRanked)$name
V(g.2)$degree <- degree(g.2)
V(gRanked)$degree <- degree(gRanked)
g.2 <- simplify(g.2)
gRanked <- simplify(gRanked)

V(g.2)$label.cex <- 1.2 * V(g.2)$degree / max(V(g.2)$degree)+ .2
V(g.2)$label.color <- rgb(0, 0, .2, .8)
V(g.2)$frame.color <- NA
egam <- (log(E(g.2)$weight)+.4) / max(log(E(g.2)$weight)+.4)
E(g.2)$color <- rgb(.5, .5, 0, egam)
E(g.2)$width <- egam
#------------------------------Clustering Analysis----------------------------------------
#Cluster Approach A
nclusters = 10
glKmeans <- kmeans(dtms, nclusters)
dis_names_dtms <- dimnames(glKmeans$centers)[[2]]

#Cluster Approach B
dt_mat <- as.matrix(dtms)
freq <- colSums(dt_mat)
ord_freq <- order(freq,decreasing = TRUE)
n_words = 100
n_words_reqd = 10
n_paragraphs =30
valid_paragraphs = which(rowSums(dt_mat[,head(ord_freq,n_words)])>n_words_reqd)
rd_mat <- dt_mat[valid_paragraphs[1:n_paragraphs],head(ord_freq,n_words)]
hc<-hclust(dist(rd_mat))
#----------------------------------Introduction Experiment------------------------------------------
Introduction<-NULL
nbrTitles<-length(as.array(Title))
Abstracts<-unlist(Abstract)
citation.Intro<-unlist(citation.Intro)
for(i in 1:nbrTitles)
{
  sentences <- as.character(Abstracts[i])
  sentences<-strsplit(sentences,split = "[\\.?!] ")
  sentences<-as.data.frame(sentences)
  nbrSentences<-nrow(sentences)
  if(nbrSentences > 3)
  {
    Introduction[i]<-str_join(citation.Intro[i],
                                           '[',as.String(i),']',
                                           'states that ',
                                           sentences[nbrSentences-2,1],
                                           sentences[nbrSentences-1,1],
                                           sentences[nbrSentences,1])
  }
}
#----------------------------------Tables--------------------------------------------------

Table.1<-xtable(DNA.df[1:10,])
Table.2<-xtable(RDAssociation.df)
#----------------------------------Figures------------------------------------------------

Figure.1<-plot(hc)
Figure.2<-plot(as.phylo(hc),cex=0.9, label.offset=1)
Figure.3<-plot(as.phylo(hc),cex=0.9, type="fan")

Figure.4<-plot(g,layout=layout.kamada.kawai, vertex.size=1, vertex.color="green")
layout.1 <- layout.fruchterman.reingold(g)
Figure.5<-plot(g.2,layout=layout.1, vertex.size=3, vertex.color="blue")
Figure.6<-plot(gRanked, layout=layout.kamada.kawai)

#---------------------------------References---------------------------------------------
Author<-xpathApply(fetch.pmc, '//PubmedArticle//Article//AuthorList/Author', function(x)
  xmlValue(xmlChildren(x)$LastName))
Authors<-xpathApply(fetch.pmc, '//PubmedArticle//Article//AuthorList', function(x)
  xmlValue(xmlChildren(x)$Author))
Title = xpathApply(fetch.pmc, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$ArticleTitle))
JournalTitle<-xpathApply(fetch.pmc, '//PubmedArticle//Article//Journal', function(x)
  xmlValue(xmlChildren(x)$Title))
JournalIssue<-xpathApply(fetch.pmc, '//PubmedArticle//Article//Journal', function(x)
  xmlValue(xmlChildren(x)$JournalIssue))

citation.BiB<-paste(Authors,":",Title,":",JournalTitle,":",JournalIssue)
#---------------------------------Function Library---------------------------------------

interactions.from.gene <- function(gene.id){                                                                                                                                                                         
  xmlrec <- entrez_fetch(db="gene", id=gene.id, rettype="xml", parsed=TRUE)                                                                                                                                          
  XML::xpathSApply(xmlrec,                                                                                                                                                                                           
                   "//Gene-commentary[Gene-commentary_heading[./text()='Interactions']]//Other-source[Other-source_src/Dbtag/Dbtag_db[./text()='GeneID']]//Other-source_anchor",
                   XML::xmlValue)                                                                                                                                  
}