library(HDMD)
library(xtable)
library(pdc)
library(lattice)
library(fracdiff)
library(tseriesEntropy)
#---------------------------------Data------------------------------------------------

grp1 <- replicate(10, arima.sim(n = 48, list(ar = c(0.8, -0.4), 
                                             ma = c(-0.22, 0.2)),
                                sd = sqrt(0.1)) )

grp2 <- replicate(10, arima.sim(n = 48, list(ar = c(-0.7, 0.1), 
                                             ma = c(0.9, -0.1)),
                                sd = sqrt(0.09)) )

#long tailed distribution
grp3 <- replicate(10, arima.sim(n = 48, list(ar = c(-0.5, 0.25), 
                                             ma = c(0.5, 0.25)),
                                rand.gen = function(n, ...) sqrt(0.1) * rt(n, df = 10 )))

#fractional 
grp4<-replicate(10,fracdiff.sim(48, ar = .2, ma = .4, d = .3))

#--------------------------------Minimum Entropic Heuristic-----------------------------
heuristic.1 <-  entropyHeuristic(grp1[,1] )
heuristic.delay.1 <-  entropyHeuristic(grp1[,1], t.min=1, t.max=6 )
heuristic.2 <-  entropyHeuristic(grp2[,1] )
heuristic.delay.2 <-  entropyHeuristic(grp2[,1], t.min=1, t.max=6 )
heuristic.3 <-  entropyHeuristic(grp3[,1] )
heuristic.delay.3 <-  entropyHeuristic(grp3[,1], t.min=1, t.max=6 )
heuristic.4 <-  entropyHeuristic(grp4[,1]$series )
heuristic.delay.4 <-  entropyHeuristic(grp4[,1]$series, t.min=1, t.max=6 )


#----------Permutation Distribution Clustering for time series---------------------------

Group.X <- cbind(grp1,grp2,grp3)

D <- pdcDist(Group.X,3)
clustering <- pdclust(Group.X)

x <- codebook(grp2[,1],m=4)
y <- codebook(grp4[,1]$series,m=4)
D.hellinger<-hellingerDistance(x,y)
D.hellinger.squared<-squaredHellingerDistance(x,y)
Divergence<-symmetricAlphaDivergence(x,y)


#----------Hierarchical cluster analysis for time series---------------------------

hc.1<-hclust(dist(Group.X))
memb.1 <- cutree(hc.1, k = 20)
memb.2 <- cutree(hc.1, k = 15)
memb.3 <- cutree(hc.1, k = 10)
cent <- NULL
for(k in 1:20){
  cent <- rbind(cent, colMeans(Group.X[memb.1 == k, , drop = FALSE]))
}
hc.2 <- hclust(dist(cent)^2, method = "cen", members = table(memb.1))
cent <- NULL
for(k in 1:15){
  cent <- rbind(cent, colMeans(Group.X[memb.2 == k, , drop = FALSE]))
}
hc.3 <- hclust(dist(cent)^2, method = "cen", members = table(memb.2))
cent <- NULL
for(k in 1:10){
  cent <- rbind(cent, colMeans(Group.X[memb.3 == k, , drop = FALSE]))
}
hc.4 <- hclust(dist(cent)^2, method = "cen", members = table(memb.3))

#----------KMeans clustering------------------------------------------------------

clustering.Kmeans<-kmeans(Group.X,3)

metric.ss <- function(x) sum(scale(x, scale = FALSE)^2)
clustering.Kmeans.fitted <- fitted(clustering.Kmeans)
head(clustering.Kmeans.fitted)
clustering.Kmeans.resid<- Group.X - fitted(clustering.Kmeans)
diagnostics.df<-cbind(clustering.Kmeans[c("betweenss", "tot.withinss", "totss")],
      c(metric.ss(clustering.Kmeans.fitted), 
        metric.ss(clustering.Kmeans.resid),    
        metric.ss(Group.X)))


#-------------------------------Tables-----------------------------------------------

Table.1<-xtable(diagnostics.df)
#------------------------------Figures----------------------------------------------
par(mfrow = c(2,1))
Figure.1<-plot(grp1[,1], type="l", 
               lty=1,col="black",
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(grp1[,2], lty=2,col="red")
lines(grp1[,3], lty=3,col="blue")
lines(grp2[,1], lty=4,col="green")
lines(grp2[,2], lty=5,col="yellow")
lines(grp2[,3], lty=6,col="orange")

Figure.2<-plot(grp3[,1], type="l", 
               lty=1,col="black",
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(grp3[,2], lty=2,col="red")
lines(grp3[,3], lty=3,col="blue")
lines(grp3[,4], lty=4,col="green")
lines(grp3[,5], lty=5,col="yellow")
lines(grp3[,6], lty=6,col="orange")

par(mfrow = c(2,2))
Figure.3A<-plot(heuristic.1)
Figure.3B<-plot(heuristic.2)
Figure.3C<-plot(heuristic.3)
Figure.3D<-plot(heuristic.4)
par(mfrow = c(2,2))
Figure.4A<-plot(heuristic.delay.1)
Figure.4B<-plot(heuristic.delay.2)
Figure.4C<-plot(heuristic.delay.3)
Figure.4D<-plot(heuristic.delay.4)

Figure.5<-levelplot(as.matrix(D), col.regions=grey.colors(20,start=0.1, end=0.95))

par(mfrow = c(1,2))
Figure.6<-plot(clustering, labels=1:30, cols=c(rep("red",10),
                                               rep("blue",10),
                                               rep("green",10)))

Figure.6A<-plot(clustering, cols=c(rep("red",10),
                                   rep("blue",10),
                                   rep("green",10)))

par(mfrow = c(2,2))
Figure.7A<-plot(hc.1)
Figure.7B<-plot(hc.2)
Figure.7C<-plot(hc.3)
Figure.7D<-plot(hc.4)

Figure.8<-plot(Group.X, col = clustering.Kmeans$cluster)
points(clustering.Kmeans$centers, col = 1:3, pch = 8, cex = 2)


#------------------------------References------------------------------------------

Reference.1<-c("",
               "",
               "")
#------------------------------Functional Library----------------------------------










