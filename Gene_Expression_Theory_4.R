library(sysBio)
library(deSolve)
library(ReacTran)
library(rootSolve)
library(fda)
library(phaseR)
library(xtable)
library(tseriesChaos)
library(corrplot)
library(plot3D)
library(scatterplot3d)
library(rgl)
library(xtable)
library(PearsonDS)
library(fitdistrplus)
#---------------------------------------------Data---------------------------------------------------------

#--------------------------------------------Model--------------------------------------------------------

Gene.Expression.Model.1 <- newModel("Gene Regulation") 
#--------------------------------------------Add reactions 
addMAreaction(Gene.Expression.Model.1, "DNA -> DNA + mRNA", r1="v1", name="Transcription")
addMAreaction(Gene.Expression.Model.1, "mRNA -> mRNA + protein", r1="v2", name="Translation")
addMAreaction(Gene.Expression.Model.1, "DNA + protein -> DNA_protein", r1="v3", name="Binding")
#-------------------------------------------Add Reaction Rates-------------------------------------------
addMAreactRate(Gene.Expression.Model.1, "v1", "assigned", "k1*DNA")
addMAreactRate(Gene.Expression.Model.1, "v2", "assigned", "k2*mRNA")
addMAreactRate(Gene.Expression.Model.1, "v3", "assigned", "k3*DNA*protein")
addMAreactRate(Gene.Expression.Model.1, "v4", "assigned", "k3r*DNA_protein")
addMAreactRate(Gene.Expression.Model.1, "v5", "assigned", "k4*mRNA")
#-------------------------------------------Add Species-------------------------------------------------
addSpecies(Gene.Expression.Model.1, "DNA", 50)
addSpecies(Gene.Expression.Model.1, "mRNA", 0)
addSpecies(Gene.Expression.Model.1, "protein", 0)
addSpecies(Gene.Expression.Model.1, "DNA_protein", 0)
#-------------------------------------------Add Parameters----------------------------------------------
addParameters(Gene.Expression.Model.1, "k1", 0.2)
addParameters(Gene.Expression.Model.1, "k2", 20)
addParameters(Gene.Expression.Model.1, "k3", 0.2)
addParameters(Gene.Expression.Model.1, "k3r", 1)
addParameters(Gene.Expression.Model.1, "k4", 1.5)
addParameters(Gene.Expression.Model.1, "k5", 1)
#-------------------------------------------Add Rules--------------------------------------------------
addRule(Gene.Expression.Model.1, "rule 1","ODEs","DNA_protein=0.5*DNA")

#-------------------------------------------Evaluate Model----------------------------------------------
makeModel(Gene.Expression.Model.1)

#------------------------------------------Simulate the Model-------------------------------------------

simResults.1 <-simulateModel(Gene.Expression.Model.1, times = seq(0, 100, by = 0.1)) 
simResults.1.stoch <- solveStoch(Gene.Expression.Model.1,10,method = "D", simName = "", tau = 0.3,
                                 f = 10, epsilon = 0.03, nc = 10) 

distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("ES_1","Protein",empMoments(simResults.1$protein)),
                                 c("ES_1 Stochastic","Protein",empMoments(simResults.1.stoch$protein)))
                                 
colnames(distribution.moments.1.df)<-c("Equation","Variable","Moment 1", "Moment 2","Moment 3","Moment 4")



#--------------------------------------------Tables-------------------------------------------------------
Table.1<-xtable(distribution.moments.1.df)


#-------------------------------------------Figures-------------------------------------------------------
Figure.1<-plotResults(simResults.1, title="Gene Expression") 
Figure.2<-plotResults(simResults.1.stoch, title="Gene Expression Stochastic")



#------------------------------------------Function Library-----------------------------------------------





