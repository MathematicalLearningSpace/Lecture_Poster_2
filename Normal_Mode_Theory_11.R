library(xtable)
library(Rknots)
library(deSolve)
#--------------------------------Data------------------------------------------

protein.5VIK<- loadProtein('5VIK')

#----------------------------Knot Analysis------------------------------------

protein.5VIK.knot <- newKnot(protein.5VIK$A)

protein.5VIK.knot.cp <- closeAndProject(protein.5VIK.knot, w = 2 )

protein.5VIK.knot.delta.p <- computeInvariant( protein.5VIK.knot.cp, invariant = 'Alexander') 
protein.5VIK.knot.jones.p <- computeInvariant( protein.5VIK.knot.cp, invariant = 'Jones', skein.sign = -1) 
protein.5VIK.knot.homfly.p <- computeInvariant( protein.5VIK.knot.cp, invariant = 'HOMFLY', skein.sign = -1) 

protein.5VIK.knot.delta.p.type<-getKnotType(polynomial = protein.5VIK.knot.delta.p, invariant = 'Alexander')
protein.5VIK.knot.jones.p.type<-getKnotType(polynomial = protein.5VIK.knot.homfly.p, invariant = 'HOMFLY')
protein.5VIK.knot.jones.p.type.full<-getKnotType(polynomial = protein.5VIK.knot.homfly.p, invariant = 'HOMFLY', full.output = TRUE)

knot.analysis.df<-data.frame()
knot.analysis.df<-rbind(c("5VIK",
                          protein.5VIK.knot.delta.p.type,
                          protein.5VIK.knot.delta.p,
                          protein.5VIK.knot.jones.p,
                          protein.5VIK.knot.homfly.p))
colnames(knot.analysis.df)<-c("Protein","Type","Alexander","Jones","HOMFLY")

#---------------------------Solution to ODE----------------------------------------
parms<-c(beta=8)
xstart <- c(x1 = 1, x2 = 1, x3 = 1)
sequence <- seq(0, 10, 0.1)
system.solution.1 <- ode(y = xstart, times = sequence,
                                  func = model.1, 
                                  parms = parms)

test<-newKnot(system.solution.1[,2:4])
test.cp<-closeAndProject(test)
c<-computeInvariant(test.cp,invariant = 'Alexander')
getKnotType(polynomial = c, invariant = 'Alexander')

test.parametric.equation<-parametric.equation(sequence)
#-------------------------------Tables-----------------------------------------

Table.1<-xtable(knot.analysis.df)

#-----------------------------Figures------------------------------------------

ramp <- colorRamp(c('blue', 'white', 'red'))
pal <- rgb( ramp(seq(0, 1, length = 109)), max = 255)

Figure.1<-plotKnot3D(protein.5VIK$A, colors = list( pal ), 
                     lwd = 8, radius = .4, showNC = TRUE, text = TRUE)
par(mfrow = c(1,2))
Figure.2<-plot(protein.5VIK.knot, main = 'original', lwd = 2)
Figure.3<-plot(protein.5VIK.cp, main = 'closed & projected', lwd = 2)
Figure.4<-plot(system.solution.1[,3], type="l", 
               lty=1,col="black", ylim=c(-50,50),
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(system.solution.1[,2], lty=2,col="black")
lines(system.solution.1[,4], lty=3,col="black")
lines(test.parametric.equation$x,lty=4,col="blue")
lines(test.parametric.equation$y,lty=5,col="red")
rug(system.solution.1[,3], side=4, col="black")
rug(system.solution.1[,2], side=4, col="black")
rug(system.solution.1[,4], side=4, col="black")

#----------------------------Linkout Resource----------------------------------

LO.1<-c("http://knotprot.cent.uw.edu.pl/view/5vik/A/")
LO.2<-c("http://katlas.org/wiki/File:Figure8knot-parametricequation.png")
#----------------------------References----------------------------------------

Reference.1<-c("Mikhail Baloban,Daria M. Shcherbakova,Sergei Pletnev,Vladimir Z. Pletnev,Clark Lagarias,Vladislav V. Verkhusha",
               "Designing brighter near-infrared fluorescent proteins: insights from structural and biochemical studies",
               "Chem. Sci.,2017,8,4546")
Reference.2<-c("ROBERT W. GHRIST and PHILIP J. HOLMES",
               "AN ODE WHOSE SOLUTIONS CONTAIN ALL KNOTS AND LINKS",
               "Int. J. Bifurcation Chaos 06, 779 (1996)")

#----------------------------Function Library----------------------------------
parametric.equation<-function(t)
{
  x=(2+cos(2*t))*cos(3*t) 
  y=(2+cos(2*t))*sin(3*t)
  return(list(x=x,y=y))
}

model.1 <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dx1 <- 7*(x2-(2/7*x1-(3/14)*(abs(x1+1)-abs(x1-1))))
    dx2 <- x1-x2+x3   
    dx3 <- -beta*x2
    res <- c(dx1, dx2, dx3)
    list(res=res)
  })
}
