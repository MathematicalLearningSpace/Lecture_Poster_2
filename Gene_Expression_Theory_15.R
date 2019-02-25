#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(protr);library(randomForest);library(pROC);library(caret)
#-----------------------------------------Data-----------------------------------------------
data(AAindex)
head(AAindex,20)

fasta.files <- list.files(patt='*.*fasta$')

extracell = readFASTA(system.file(
  "protseq/extracell.fasta", package = "protr"))
mitonchon = readFASTA(system.file(
  "protseq/mitochondrion.fasta", package = "protr"))

AT.mitonchon=readFASTA("ArabidopsisThalianaMitochondrion.fasta")

extracell = extracell[(sapply(extracell, protcheck))]
mitonchon = mitonchon[(sapply(mitonchon, protcheck))]

labels = as.factor(c(rep(0, 
                         length(extracell)), 
                     rep(1, 
                         length(mitonchon))))
#----------------------------------------Data Set 1------------------------------------------
x1 = t(sapply(extracell, extractAPAAC))
x2 = t(sapply(mitonchon, extractAPAAC))
#-----------Compute Composition, Transition, and Distribution Descriptors--------------------
x3 = t(sapply(extracell,extractCTDC))
x4 = t(sapply(extracell,extractCTDT))
x5 = t(sapply(extracell,extractCTDD))

X.1 = rbind(x1, x2)

X.1.train.idx = c(
  sample(1:nrow(x1), round(nrow(x1) * 0.75)),
  sample(nrow(x1) + 1:nrow(x2), round(nrow(x2) * 0.75))
)
X.1.test.idx = setdiff(1:nrow(X.1), X.1.train.idx)
X.1.train = X.1[X.1.train.idx, ]
X.1.test = X.1[X.1.test.idx, ]

y.1.train = labels[X.1.train.idx]
y.1.test = labels[X.1.test.idx]
#----------------------------------------Data Set 2-----------------------------------------
prop.1<-runif(20,0,10)
properties.modify = data.frame(
  AccNo = c("Prop_1", "Prop_2", "Prop_3","Prop_4"),
  A = c(0.62,  -0.5, 15,prop.1[1]),  R = c(-2.53,   3, 101,prop.1[11]),
  N = c(-0.78,  0.2, 58,prop.1[2]),  D = c(-0.9,    3, 59,prop.1[12]),
  C = c(0.29,    -1, 47,prop.1[3]),  E = c(-0.74,   3, 73,prop.1[13]),
  Q = c(-0.85,  0.2, 72,prop.1[4]),  G = c(0.48,    0, 1,prop.1[14]),
  H = c(-0.4,  -0.5, 82,prop.1[5]),  I = c(1.38, -1.8, 57,prop.1[15]),
  L = c(1.06,  -1.8, 57,prop.1[6]),  K = c(-1.5,    3, 73,prop.1[16]),
  M = c(0.64,  -1.3, 75,prop.1[7]),  F = c(1.19, -2.5, 91,prop.1[17]),
  P = c(0.12,     0, 42,prop.1[8]),  S = c(-0.18, 0.3, 31,prop.1[18]),
  T = c(-0.05, -0.4, 45,prop.1[9]),  W = c(0.81, -3.4, 130,prop.1[19]),
  Y = c(0.26,  -2.3, 107,prop.1[10]), V = c(1.08, -1.5, 43,prop.1[20]))

APAAC.modify<-function(x) extractAPAAC(x, props = c("Hydrophobicity", "Hydrophilicity",
                                                "CIDH920105", "BHAR880101",
                                                "CHAM820101", "CHAM820102",
                                                "Prop_1", "Prop_2", "Prop_3","Prop_4"),
             customprops = properties.modify)

x6 = t(sapply(extracell[1:20], APAAC.modify))
x7 = t(sapply(mitonchon[1:20], APAAC.modify))

X.2  = rbind(x6, x7)

X.2.train.idx = c(
  sample(1:nrow(x6), round(nrow(x6) * 0.75)),
  sample(nrow(x7) + 1:nrow(x7), round(nrow(x7) * 0.75))
)
X.2.test.idx = setdiff(1:nrow(X.2), X.2.train.idx)
X.2.train = X.2[X.2.train.idx, ]
X.2.test = X.2[X.2.test.idx, ]

labels.2 = as.factor(c(rep(0, 
                         length(extracell[1:20])), 
                     rep(1, 
                         length(mitonchon[1:20]))))

y.2.train = labels.2[X.2.train.idx]
y.2.test = labels.2[X.2.test.idx]

#---------------------------Machine Learning Classification Model--------------------------------------------------

rf.fit.1 = randomForest(X.1.train, y.1.train, cv.fold = 5)
rf.predict.1 = predict(rf.fit.1, newdata = X.1.test, type = "prob")[, 1]

rf.fit.2 = randomForest(X.2.train, y.2.train, cv.fold = 5)
rf.predict.2 = predict(rf.fit.2, newdata = X.2.test, type = "prob")[, 1]

#-----------------------------------------Tables----------------------------------------------

#Table.1<-xtable()

#----------------------------------------Figures----------------------------------------------

Figure.1<-plot.roc(y.1.test, rf.predict.1, grid = TRUE, print.auc = TRUE)
Figure.2<-plot.roc(y.2.test, rf.predict.2, grid = TRUE, print.auc = TRUE)

#----------------------------------------References------------------------------------------

Reference.1<-c("Kuo-Chen Chou.", 
"Prediction of Protein Cellular Attributes Using Pseudo-Amino Acid Composition.", 
"PROTEINS: Structure, Function, and Genetics, 2001, 43: 246-255.")

#----------------------------------------------Function Library------------------------------------
#--------Function Template Library for Classroom Presentation and Modification---------------------
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
