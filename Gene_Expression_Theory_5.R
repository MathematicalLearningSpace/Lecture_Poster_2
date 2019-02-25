#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------------------R API --------------------------------------------
library(fractal);library(fractaldim);library(fracdiff);library(timeSeries);library(RandomFields);library(fNonlinear)
library(stats);library(boot);library(vars);library(urca);library(xtable);library(gridExtra);library(wmtsa);library(stringr)
library(stringi);library(Hmisc);library(Kendall);library(Deducer);library(corrplot);library(PearsonDS)
library(dplyr);library(readr);library(readxl)

#---------------------------------------Data Sets for Classroom----------------------------------------
#----------------------------Formatted Student Notes for the Classroom Lecture-----------
article.Notes<-readlines('Notes_Gene_Expression_Theory_5.txt')
AvailableResources <- read_csv("AvailableResources.txt")
#----------------------------Filter and Sort the Features------------------------------------------------
View(AvailableResources)
#---------------------------------------Gene By Cell Cycle Phase------------------------------------
#----------------------------Formatted Student Notes for the Classroom Lecture-----------
Gene.Cell.Cycle <- read_delim("Gene_Cell_Cycle.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
View(Gene.Cell.Cycle)

Gene.Cell.Cycle.S_Phase<-Gene.Cell.Cycle$`Gene Symbol`[grep('S phase',Gene.Cell.Cycle$`Assigned Phase`)]

Gene.Cell.Cycle.Phase.df<-data.frame()
Gene.Cell.Cycle.Phase.df<-cbind(c(Gene.Cell.Cycle.S_Phase))
colnames(Gene.Cell.Cycle.Phase.df)<-c("S-Phase")

#---------------------------------------Gene Expression Data Set---------------------------------------------
GeneExpression <- read.delim("dataPlusScores_all5.txt")
Gene.Description<-as.array(GeneExpression[,2])
#-----------RNA collected for points (typically every 1-2 h) for 30 h (Thy-Thy1),------------------ 
#-----------44 h (Thy-Thy2), 46 h (Thy-Thy3), 36 h (Thy-Noc), or 14 h (Shake) '--------------------
ThyThy1<-c(13:25)
gene.P53.ThyThy1<-t(GeneExpression[8327,13:25])
ThyThy2<-c(26:51)
gene.P53.ThyThy2<-t(GeneExpression[8327,26:51])
ThyThy3<-c(53:100)
gene.P53.ThyThy3<-t(GeneExpression[8327,53:100])
ThyNOC<-c(101:120)
gene.P53.ThyNOC<-t(GeneExpression[8327,101:120])
Shake<-c(121:130)
gene.P53.Shake<-t(GeneExpression[8327,121:130])
#------------- (a) phase (0-4 h)- (b) entered G2 phase (5-6 h)-synchronous mitosis at (c) 7-8 h,then--- 
#-------------(d) S phase, then (e) full cell cycle at 14-16 h -------------------------------------
#-------------(f) three cell cycles for 48 hours---------------------------------------------------- 
P53.Description<-Gene.Description[grep('P53',Gene.Description)]
P53<-t(GeneExpression[8327,53:100])
P53.1<-t(GeneExpression[3291,53:100])
colnames(P53)<-c('TP53')
MDM2.Description<-Gene.Description[grep('MDM2',Gene.Description)]
MDM2<-t(GeneExpression[2610,53:100])
WDR5.Description<-grep('WDR5',Gene.Description)
WDR5<-t(GeneExpression[9545,53:100])
#--------------------------------------Gene Study---------------------------------------------------

Gene.Study<-as.matrix(cbind(P53,P53.1,MDM2,WDR5))
Gene.Names<-c('P53','P53.1','MDM2','WDR5')

Gene.Study.df<-as.data.frame(Gene.Study)
names(Gene.Study.df)<-Gene.Names

#--------------------------------------Summary Statistics--------------------------------------------
Gene_Stats.df<-data.frame()
Gene_Stats.df<-rbind(c(gene.statistics.descriptive(P53)),
                          c(gene.statistics.descriptive(P53.1)),
                          c(gene.statistics.descriptive(MDM2)),
                          c(gene.statistics.descriptive(WDR5)))
colnames(Gene_Stats.df)<-c("Min","Max","Mean","Std","Skew","Kurtosis")
rownames(Gene_Stats.df)<-c(Gene.Description[8327],
                           Gene.Description[3291],
                           Gene.Description[2610],
                           Gene.Description[9545])

#-------------------------------------Distribution Statistics---------------------------------------
Gene_Dist.df<-data.frame()
Gene_Dist.df<-rbind(c(gene.statistics.distribution(P53)),
                          c(gene.statistics.distribution(P53.1)),
                          c(gene.statistics.distribution(MDM2)),
                          c(gene.statistics.distribution(WDR5)))
colnames(Gene_Dist.df)<-c("Type","Location","Scale")
rownames(Gene_Dist.df)<-c(Gene.Description[8327],
                          Gene.Description[3291],
                          Gene.Description[2610],
                          Gene.Description[9545])

#--------------------------------------Gene Correlation-----------------------------------------------
gene.correlation <- cor(Gene.Study)

#---------------------------------------Tables-------------------------------------------------------

Table.1<-xtable(Gene_Stats.df)
Table.2<-xtable(Gene_Dist.df)

#---------------------Figures for Classroom Discussion-------------------------------------------------------

Figure.1<-plot(P53,type = 'l', col = "black", lwd = 2, 
               lty=1,main = "TP53, MDM2 and WDR5 Expression Cell Cycle",ylab="Expression",xlab="Hours")
lines(P53.1,lty=2,col="red")
lines(MDM2,lty=3,col="blue")
lines(WDR5,lty=3,col="green")
grid()
legend("bottomright", col = c("black","red","blue","green"), 
       lty = 1:1, cex=0.75,
       legend = c("p53","P53.1","MDM2","WDR5"))

Figure.2<-corrplot(gene.correlation)
title(main="Gene Correlation Plot")
#-------Function Library to be Modfied by Students in the Classroom------------------------------------------------

gene.statistics.descriptive <-function(x)
{
  stat.df<-NULL
  stat.df.2<-NULL
  for(i in 1:ncol(x))
  {
    x.min<-min(x[,i])
    x.max<-max(x[,i])
    x.mean<-mean(x[,i])
    x.std<-var(x[,i])^(1/2)
    x.skew<-skewness(x[,i])
    x.kurtosis<-kurtosis(x[,i])
    stat.df<-rbind(x.min,
                   x.max,
                   x.mean,
                   x.std,
                   x.skew,
                   x.kurtosis)
    stat.df.2<-cbind(stat.df.2,stat.df)
    stat.df<-NULL
  }
  return(stat.df.2)
}
gene.statistics.distribution<-function(x){
  stat.df<-NULL
  stat.df.2<-NULL
  for(i in 1:ncol(x))
  {
    x.mean<-mean(x[,i])
    x.std<-var(x[,i])^(1/2)
    x.skew<-skewness(x[,i])
    x.kurtosis<-kurtosis(x[,i])
    moments <- c(mean=x.mean,variance=x.std^2,skewness=x.skew,kurtosis=x.kurtosis)
    #distribution Test
    fit<-pearsonFitML(x[,1])
    #fit.mm<-pearsonFitM(moments=moments)
    stat.df<-rbind(fit$type,
                   fit$location,
                   fit$scale)
                   #fit.mm$type)
    
    stat.df.2<-cbind(stat.df.2,stat.df)
    stat.df<-NULL
  }
  return(stat.df.2)
}
