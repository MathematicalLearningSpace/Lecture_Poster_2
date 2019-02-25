#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(xtable);library(qdap);library(qdapDictionaries);library(qdapRegex);library(qdapTools);library(stringi);library(text2vec);library(readability)
library(syllable);library(httr);library(XML);library(psych);library(Hmisc)

#-----------------------------------------Data Set Tables--------------------------------------------------
url <- "https://en.wikipedia.org/wiki/Ligand"
url.1<-"https://en.wikipedia.org/wiki/Statistical_data_type"
url.2<-"https://en.wikipedia.org/wiki/Graph_theory"
url.3<-"https://en.wikipedia.org/wiki/Data_visualization"
url.4<-"https://en.wikipedia.org/wiki/Mitochondrial_DNA"
url.5<-"https://en.wikipedia.org/wiki/Harmonic_oscillator"

url.content <- GET(url.4)
document <- readHTMLTable(doc=content(url.content, "text"))
document[1]

document.mitochondrial.DNA<-as.data.frame(document[1])
write.csv(document.mitochondrial.DNA, file="Mitochondrial_DNA.csv",row.names = F)

document.ligand<-as.data.frame(document[2])
write.csv(document.ligand, file="ligand.csv",row.names = F)

Table.1.data.set <- data.frame(y=rnorm(5), 
                               x1=c(1:5), 
                               x2=c(TRUE, TRUE, FALSE, FALSE, FALSE),
                               X3=letters[1:5])

Table.1.data.set.description <- c("The value of y ranges from the minimum value to the maximum value. ",
                                  "The value of x1 ranges from the minimum value to the maximum value. ",
                                  "The value of x2 has n true values and (1-n) false values. ",
                                  "The value of x3 has k categories. ",
                                  "This data set is from this data reference 1. "
)

Table.1.statistics<-stri_stats_latex(Table.1.data.set.description)

Table.1.paragraph<-stri_join(Table.1.data.set.description[1],
                             Table.1.data.set.description[2],
                             Table.1.data.set.description[3],
                             Table.1.data.set.description[4],
                             Table.1.data.set.description[5])

Table.1.paragraph.token.1<-itoken(Table.1.paragraph, tokenizer = word_tokenizer,
                                  progressbar = T)

vocab.Table.1<-create_vocabulary(Table.1.paragraph.token.1, ngram = c(1L, 1L)) 
vocab.Table.2<-create_vocabulary(Table.1.paragraph.token.1, ngram = c(2L, 2L))
vocab.Table.3<-create_vocabulary(Table.1.paragraph.token.1, ngram = c(3L, 3L))

Figure.1.data.set.description<-c("Figure 1 is line graph for data frame name .",
                                 "shows the relationship between Y on the X axis among .",
                                 "The value of X2 .",
                                 "The value of X3 .",
                                 "The value of X4 .",
                                 "It shows the frequency and scale of changes in the first and second partial derivatives. ")

Figure.1.statistics<-stri_stats_latex(Figure.1.data.set.description)

Figure.1.paragraph<-stri_join(Figure.1.data.set.description[1],
                              Figure.1.data.set.description[2],
                              Figure.1.data.set.description[3],
                              Figure.1.data.set.description[4],
                              Figure.1.data.set.description[5])

Figure.1.paragraph.token.1<-itoken(Figure.1.paragraph, tokenizer = word_tokenizer,
                                   progressbar = T)

vocab.Figure.1<-create_vocabulary(Figure.1.paragraph.token.1, ngram = c(1L, 1L))
vocab.Figure.2<-create_vocabulary(Figure.1.paragraph.token.1, ngram = c(2L, 2L))
vocab.Figure.3<-create_vocabulary(Figure.1.paragraph.token.1, ngram = c(3L, 3L))

#---------------------------------Analysis of the Sentences---------------------------------------

Table.readability.scores.df<-readability(Table.1.paragraph,Table.1.paragraph)
Figure.readability.scores.df<-readability(Figure.1.paragraph,Figure.1.paragraph)


x<-rnorm(100)
print(moment.1(x, "mean"))
print(moment.2(x, "median"))
print(moment.3(x, "trimmed"))


#-------------------------------References--------------------------------------------------------

Reference.1<-stringr::str_c("Wikipedia contributors,"
               "Harmonic oscillator.",
               "Wikipedia, The Free Encyclopedia, 25 Feb. 2019.")


#-----------Function Library To be Designed and Testing in the Classroom---------------------------------------------------

moment.1 <- function(x, measure) {
  switch(measure,
         mean = {
           total.x<-sum(x)
           mean = total.x/length(x)
         },
         median = median(x),
         trimmed = mean(x, trim = .15)
  )
}

Figure.choice.template<-function(x,figureChoice)
{
  switch(figureChoice,
         A={
           #Line Graph
         },
         B={
           #Scatterplot matrix
         },
         C={
           #histogram
         },
         D={
           #Graph 
         })
  
}

Table.choice.template<-function(x,tableChoice)
{
  switch(tableChoice,
         A={
           #Data Frame
         },
         B={
           #Matrix
         },
         C={
           #Summary Statistics
         },
         D={
           #Generic Table
         })
  
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
