library(xtable)
library(Rknots)
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


#----------------------------Linkout Resource----------------------------------

LO.1<-c("http://knotprot.cent.uw.edu.pl/view/5vik/A/")
LO.2<-c("http://katlas.org/wiki/File:Figure8knot-parametricequation.png")
#----------------------------References----------------------------------------

Reference.1<-c("Mikhail Baloban,Daria M. Shcherbakova,Sergei Pletnev,Vladimir Z. Pletnev,Clark Lagarias,Vladislav V. Verkhusha",
               "Designing brighter near-infrared fluorescent proteins: insights from structural and biochemical studies",
               "Chem. Sci.,2017,8,4546")


#----------------------------Function Library----------------------------------
