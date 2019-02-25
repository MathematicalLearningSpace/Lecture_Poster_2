#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(xtable);library(text2vec);library(data.table);library(magrittr);library(glmnet);library(rentrez);library(readr);library(LDAvis)

#------------------------------------------Data-----------------------------------------------
article.files <- list.files(patt='*.*csv$')
#-------------------------------Formatted Student Notes------------------------------------------
print(article.files)
#[1] "publications_AML.csv"             
#[2] "publications_Apoptosis.csv"       
#[3] "publications_AR.csv"              
#[4] "publications_BladderCancer.csv"   
#[5] "publications_BreastCancer.csv"    
#[6] "publications_Calcium.csv"         
#[7] "publications_Cell_Cycle_G1.csv"   
#[8] "publications_CML.csv"             
#[9] "publications_Colorectal.csv"      
#[10] "publications_Endometrial.csv"     
#[11] "publications_GastricCancer.csv"   
#[12] "publications_Hepatocellular.csv"  
#[13] "publications_HH.csv"              
#[14] "publications_HIF_1.csv"           
#[15] "publications_JAK_STAT.csv"        
#[16] "publications_KEAP1_NRF2.csv"      
#[17] "publications_MAPK_ERK.csv"        
#[18] "publications_NonsmallCellLung.csv"
#[19] "publications_NOTCH.csv"           
#[20] "publications_Other_RAS.csv"       
#[21] "publications_Pancreatic.csv"      
#[22] "publications_PI3K.csv"            
#[23] "publications_Prostate.csv"        
#[24] "publications_smallCellLung.csv"   
#[25] "publications_Telomerase.csv"      
#[26] "publications_TGFB.csv"            
#[27] "publications_ThyroidCancer.csv"   
#[28] "publications_WNT.csv"             
z.df<-data.frame()
z.df<-read_csv(article.files[25])
x<-rbinom(nrow(z.df),1,0.5)
z.df<-cbind(z.df,score=x)
#-----------------------Filter, Search and Sort Data Set---------------------------
View(z.df)

Data_Dictionary <- read_csv("Data_Dictionary.txt")

#-----------------------------------------Process Titles Data----------------------------------
Titles<-z.df$Title

#-----------------------------------------Process Abstracts Data-------------------------------
z.dt<-data.table(z.df)
#setDT(z.dt)
setkey(z.dt,"Pubmed ID")
set.seed(2017L)
pubmed.all_ids = z.df$`Pubmed ID`

train.ids = sample(pubmed.all_ids,0.75*nrow(z.df))
test.ids = setdiff(pubmed.all_ids, train.ids)
train = z.dt[J(train.ids)]
test = z.dt[J(test.ids)]

prep.fun = tolower
tok.fun = word_tokenizer

Abstract.train.tokens = train$Abstract %>% 
  prep.fun %>% 
  tok.fun

Title.train.tokens = train$Title %>% 
  prep.fun %>% 
  tok.fun

Abstract.it.train = itoken(Abstract.train.tokens, 
                  ids = train$`Pubmed ID`,
                  progressbar = FALSE)

Abstract.it.test = test$Abstract %>% 
  prep.fun %>% tok.fun %>% 
  itoken(ids = test$`Pubmed ID`, progressbar = FALSE)

Title.it.train = itoken(Title.train.tokens, 
                           ids = train$`Pubmed ID`,
                           progressbar = FALSE)

Title.it.test = test$Title %>% 
  prep.fun %>% tok.fun %>% 
  itoken(ids = test$`Pubmed ID`, progressbar = FALSE)


#--------------------------Vocabulary Work on Abstracts and Titles-----------------------------

vocab = create_vocabulary(Abstract.it.train)
vocab
Title.vocab=create_vocabulary(Title.it.train)
Title.vocab

stop_words = c("i", "me", "my", "myself", 
               "we", "our", "ours", "ourselves", 
               "you", "your", "yours")

vocab = create_vocabulary(Abstract.it.train, 
                          stopwords = stop_words)
vocab
vectorizer = vocab_vectorizer(vocab)

#-------------------------Prune the Absrtract Vocabulary Tree-------------------------------------------

pruned.vocab = prune_vocabulary(vocab, 
                                term_count_min = 10, 
                                doc_proportion_max = 0.5,
                                doc_proportion_min = 0.001)

vectorizer.pruned = vocab_vectorizer(pruned.vocab)

#------------Tokenize the Title and Abstract Vocabulary---------------------------------------------

token.Title<-word_tokenizer(tolower(z.dt$Title))
vocab.Title<-create_vocabulary(itoken(token.Title),ngram=c(1L,3L))

vectorizer.Title<-vocab_vectorizer(vocab.Title)

tokens = word_tokenizer(tolower(z.dt$Abstract))
v = create_vocabulary(itoken(tokens))
v = prune_vocabulary(v, term_count_min = 3, 
                     doc_proportion_max = 0.5)

it = itoken(tokens)
vectorizer.tokens = vocab_vectorizer(v)

vocab.grams.2 = create_vocabulary(Abstract.it.train, 
                                  ngram = c(1L, 2L))
vocab.grams.2 = prune_vocabulary(vocab.grams.2, 
                                 term_count_min = 3, 
                         doc_proportion_max = 0.75)
vocab.grams.2

vocab.grams.3 = create_vocabulary(Abstract.it.train, 
                                  ngram = c(1L, 3L))
vocab.grams.3 = prune_vocabulary(vocab.grams.3, 
                                 term_count_min = 3, 
                                 doc_proportion_max = 0.75)
vocab.grams.3

#------------------------Bi and Tr-Gram Approach to the Vocabulary--------------------------

bigram_vectorizer = vocab_vectorizer(vocab.grams.2)
trigram_vectorizer=vocab_vectorizer(vocab.grams.3)

#-----------------------------------------Models----------------------------------------------
dtm_train = create_dtm(Abstract.it.train, vectorizer)
dtm_test = create_dtm(Abstract.it.test, vectorizer)

tfidf = TfIdf$new()
dtm_train_tfidf = fit_transform(dtm_train, tfidf)
dtm_test_tfidf = create_dtm(Abstract.it.test, vectorizer)
dtm_test_tfidf = transform(dtm_test_tfidf, tfidf)

NFOLDS = 4
glmnet_classifier.1 = cv.glmnet(x = dtm_train_tfidf, 
                                y = train[['score']], 
                              family = 'binomial', 
                              alpha = 1,
                              type.measure = "auc",
                              nfolds = NFOLDS,
                              thresh = 1e-3,
                              maxit = 1e3)

print(paste("max AUC =", round(max(glmnet_classifier.1$cvm), 4)))
preds.1 = predict(glmnet_classifier.1, dtm_test_tfidf, type = 'response')[,1]
glmnet:::auc(test$score, preds.1)

glmnet_classifier.2 = cv.glmnet(x = dtm_train, y = train[['score']], 
                              family = 'binomial', 
                              alpha = 1,
                              type.measure = "auc",
                              nfolds = NFOLDS,
                              thresh = 1e-3,
                              maxit = 1e3)
print(paste("max AUC =", round(max(glmnet_classifier.2$cvm), 4)))
preds.2 = predict(glmnet_classifier.2, dtm_test, type = 'response')[,1]
glmnet:::auc(test$score, preds.2)

#-------------------------Vocabulary Changes I------------------------

vectorizer.pruned = vocab_vectorizer(pruned.vocab)

dtm_train  = create_dtm(Abstract.it.train, vectorizer.pruned)
glmnet_classifier.3 = cv.glmnet(x = dtm_train, y = train[['score']], 
                              family = 'binomial', 
                              alpha = 1,
                              type.measure = "auc",
                              nfolds = NFOLDS,
                              thresh = 1e-3,
                              maxit = 1e3)
print(paste("max AUC =", round(max(glmnet_classifier.3$cvm), 4)))
dtm_test = create_dtm(Abstract.it.test,vectorizer.pruned)
preds.3 = predict(glmnet_classifier.3, dtm_test, type = 'response')[,1]
glmnet:::auc(test$score, preds.3)

#------------------------Vocabulary Changes II-------------------------

dtm_train = create_dtm(Abstract.it.train, bigram_vectorizer)
glmnet_classifier.4 = cv.glmnet(x = dtm_train, y = train[['score']], 
                              family = 'binomial', 
                              alpha = 1,
                              type.measure = "auc",
                              nfolds = NFOLDS,
                              thresh = 1e-3,
                              maxit = 1e3)
print(paste("max AUC =", round(max(glmnet_classifier.4$cvm), 4)))

dtm_test = create_dtm(Abstract.it.test, bigram_vectorizer)
preds.4 = predict(glmnet_classifier.4, dtm_test, type = 'response')[,1]
glmnet:::auc(test$score, preds.4)

#-------------------------------Relaxed Word Movers Model-------------------------------

dtm = create_dtm(Abstract.it.train, vectorizer.tokens)
tcm = create_tcm(Abstract.it.train, vectorizer.tokens, skip_grams_window = 5)

glove_model = GloVe$new(word_vectors_size = 50, vocabulary = v, x_max = 10)
wv = glove_model$fit_transform(tcm, n_iter = 10)
wv = wv + t(glove_model$components)

rwmd_model = RWMD$new(wv)
rwmd_dist = dist2(dtm[1:100, ], dtm[1:10, ], 
                  method = rwmd_model, norm = 'none')

#-----------------------------Binormal separation----------------------------------------

model_bns = BNS$new()
dtm_bns = model_bns$fit_transform(dtm, head(z.dt$score, 100))

#----------------------------LDA Model--------------------------------------------------
n.topics<-10
dtm = create_dtm(Abstract.it.train, vectorizer.tokens)

lda_model = LDA$new(n.topics)

doc_topic_distr = lda_model$fit_transform(dtm, n_iter = 20)
doc_topic_word_distr_10<-lda_model$topic_word_distribution

lda_model$topic_word_distribution
#lda_model$components

lda_model.Top.10<-lda_model$get_top_words(n=10)

lda_model.perplexity<-perplexity(dtm,doc_topic_word_distr_10,doc_topic_distr)

#----------------------------LSA Model---------------------------------------------------

lsa_model = LatentSemanticAnalysis$new(n.topics,method="randomized")

doc_topic_distr.1 = lsa_model$fit_transform(dtm)
doc_topic_distr.2 = fit_transform(dtm, lsa_model)
lsa_model$components

#-------------------------------Model Comparisions---------------------------------------
Model.Comparision.df<-data.frame()
Model.Comparision.df<-rbind(c("1",round(max(glmnet_classifier.1$cvm), 4),glmnet:::auc(test$score, preds.1)),
                            c("2",round(max(glmnet_classifier.2$cvm), 4),glmnet:::auc(test$score, preds.2)),
                            c("3",round(max(glmnet_classifier.3$cvm), 4),glmnet:::auc(test$score, preds.3)),
                            c("4",round(max(glmnet_classifier.4$cvm), 4),glmnet:::auc(test$score, preds.4))
)
colnames(Model.Comparision.df)<-c("Model","AUC max","AUC")

#------------------------------LDA and LSA Analysis----------------------------------------

Semantic.Analysis.df<-data.frame()
Semantic.Analysis.df<-rbind(c(lda_model.perplexity)
)
colnames(Semantic.Analysis.df)<-c("Perplexity")

#-----------------------------Topic Key word Tables for PubMed searches--------------------

Topic.Key.Word.df<-data.frame()
Topic.Key.Word.df<-cbind(c(lda_model.Top.10[,1]),c(lda_model.Top.10[,2]),c(lda_model.Top.10[,3]))
colnames(Topic.Key.Word.df)<-c("LDA_Top_10 A","LDA_Top_10 B","LDA_Top_10 C")


#----------------------------Title Analysis------------------------------------------------

Title.Analysis.df<-data.frame()
Title.Analysis.df<-cbind(z.df$Title,z.df$score)
colnames(Title.Analysis.df)<-c("Title","Score")

#-----------------------------------------Tables---------------------------------------------

Table.1<-xtable(Model.Comparision.df)
Table.1
Table.2<-xtable(Semantic.Analysis.df)
Table.2
Table.3<-xtable(Topic.Key.Word.df)
Table.3
Table.4<-xtable(Title.Analysis.df)
Table.4
#----------------------------------------Figures---------------------------------------------
par(mfcol = c(2, 2))
Figure.1<-plot(glmnet_classifier.1,main="Model 1")
Figure.2<-plot(glmnet_classifier.2,main="Model 2")
Figure.3<-plot(glmnet_classifier.3,main="Model 3")
Figure.4<-plot(glmnet_classifier.4,main="Model 4")

par(mfcol = c(1, 1))
Figure.5<-plot(rwmd_dist[1,],type="l",main="RWMD")
for(i in 2:nrow(rwmd_dist))
{
  lines(rwmd_dist[i,],lty=i)
}
Figure.6<-lda_model$plot()

#----------------------------------------References-----------------------------------------

Reference.1<-c("Jonathan Chang (2012).",
               "lda: Collapsed Gibbs sampling methods for topic models.", 
               "R package version 1.3.2. http://CRAN.R-project.org/package=lda.")

Reference.2<-c("Blei, David M., Ng, Andrew Y., and Jordan, Michael I. (2003).", 
                "Latent Dirichlet Allocation",
                "Journal of Machine Learning Research, Volume 3, pages 993-1022.")

Reference.3<-c("Sievert, C. and Shirley, K. (2014)",
               "LDAvis: A Method for Visualizing and Interpreting Topics",
               "ACL Workshop on Interactive Language Learning, Visualization, and Interfaces. http://nlp.stanford.edu/events/illvi2014/papers/sievert-illvi2014.pdf")


#----------------------------------------Function Library-----------------------------------
