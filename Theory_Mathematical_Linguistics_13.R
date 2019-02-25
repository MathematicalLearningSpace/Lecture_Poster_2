#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(rentrez);library(xtable);library(stringi);library(readr)
#------------------------------------------Data-------------------------------------------------------------------------

KEGG.IDs.df <- as.data.frame(read_csv("data/KEGG_IDs_For_Query.txt"))
KEGG.Cancer.IDs.df<- as.data.frame(read_csv("data/KEGG_Cancer_IDs_For_Query.txt"))

View(KEGG.IDs.df)

#------------------------------------------KEGG Data Processing---------------------------------------------------------
KEGG.Databases<-listDatabases()
KEGG.organisms<-keggList("organism")
KEGG.Brite<-keggList("brite")

#------------------------------------------Genetic information processing-----------------------------------------------
KEGG.organisms.df<-as.data.frame(KEGG.organisms)
KEGG.organisms.df$species
KEGG.organisms.df$species[210]

KEGG.DNA.repair<-keggGet("br:ko03029")
KEGG.DNA.repair.1<-keggGet("path:pxb03420")
KEGG.DNA.repair.2<-keggGet("path:pxb03440")

KEGG.Proteasome<-keggGet("br:ko03051")

KEGG.Proteasome.20s<-keggGet("path:pxb03050")

KEGG.Ubiquitin.System<-keggGet("br:ko04121")
KEGG.Ubiquitin.System.1<-keggGet("path:pxb04120")
KEGG.Ubiquitin.System.2<-keggGet("path:pxb03420")

KEGG.organisms.AT<-dfKEGG.organisms.df[177,]

KEGG.pathway.Database<-keggList(KEGG.Databases[1])
KEGG.pathways<-keggLink("pathway", "hsa")
KEGG.pathways.List<-as.list(keggList("pathway", "hsa"))
KEGG.pathways.04010<-keggGet("path:hsa04010")
KEGG.pathways.04011<-keggGet("path:map04011")

KEGG.pathways.cancer.basal.cell.carcinoma<-keggGet("path:hsa05217")
KEGG.pathways.cancer.pathways<-keggGet("path:hsa05200")

KEGG.IDs.signaling.df<-KEGG.IDs.df[grep("signaling",KEGG.IDs.df$Description),]
KEGG.IDs.signaling.Wnt.df<-KEGG.IDs.df[grep("Wnt",KEGG.IDs.df$Description),]
KEGG.pathways.04310<-keggGet("path:map04310")

#----------------------------------------KEGG Signaling Pathways : Gene, Articles---------------------------------------
number.even.1<-seq(0,length(KEGG.pathways.04010[[1]]$GENE),2)
KEGG.pathways.04010.Gene<-KEGG.pathways.04010[[1]]$GENE[number.even.1]
KEGG.pathways.04010.Articles<-KEGG.pathways.04010[[1]]$REFERENCE

KEGG.pathways.04010.Articles.References<-NULL
for(i in 1:length(KEGG.pathways.04010.Articles))
{
KEGG.pathways.04010.Articles.References[i]<-stri_join(KEGG.pathways.04010.Articles[[i]]$REFERENCE,";",
          KEGG.pathways.04010.Articles[[i]]$AUTHORS,";",
          KEGG.pathways.04010.Articles[[i]]$TITLE,";",
          KEGG.pathways.04010.Articles[[i]]$JOURNAL)
}
KEGG.pathways.04010.Articles.References.df<-as.data.frame(KEGG.pathways.04010.Articles.References)
colnames(KEGG.pathways.04010.Articles.References.df)<-c("Reference")

#--------------------------------------KEGG Cancer Pathways : Gene, Articles-----------------------------------------------
number.even.2<-seq(0,length(KEGG.pathways.cancer.pathways[[1]]$GENE),2)
KEGG.pathways.cancer.pathway.Gene.Names<-KEGG.pathways.cancer.pathways[[1]]$GENE[number.even.2]

KEGG.pathways.cancer.pathway.Gene.proteins.heatshock<-KEGG.pathways.cancer.pathway.Gene.Names[grep("heat shock protein",KEGG.pathways.cancer.pathway.Gene.Names)]
search.term<-stri_extract_all_words(KEGG.pathways.cancer.pathway.Gene.proteins.heatshock[1])

publications.proteins.heatshock <- entrez_search(db="pubmed", term=search.term[[1]][1], retmax=40)

KEGG.pathways.cancer.pathway.Articles<-KEGG.pathways.cancer.pathways[[1]]$REFERENCE
KEGG.pathways.cancer.pathway.Articles.References<-NULL
for(i in 1:length(KEGG.pathways.cancer.pathway.Articles))
{
  KEGG.pathways.cancer.pathway.Articles.References[i]<-stri_join(KEGG.pathways.cancer.pathway.Articles[[i]]$REFERENCE,";",
                                                        KEGG.pathways.cancer.pathway.Articles[[i]]$AUTHORS,";",
                                                        KEGG.pathways.cancer.pathway.Articles[[i]]$TITLE,";",
                                                        KEGG.pathways.cancer.pathway.Articles[[i]]$JOURNAL)
}
KEGG.pathways.cancer.pathway.Articles.References.df<-as.data.frame(KEGG.pathways.cancer.pathway.Articles.References)
colnames(KEGG.pathways.cancer.pathway.Articles.References.df)<-c("Reference")

#-----------------------------------------Tables-------------------------------------------

Table.1<-xtable(KEGG.IDs.df)
Table.2<-xtable(KEGG.Cancer.IDs.df)
Table.3<-xtable(KEGG.IDs.signaling.df)
Table.4<-xtable(KEGG.pathways.04010.Articles.References.df)
Table.5<-xtable(head(KEGG.pathways.04010.Articles.References.df))

#-----------------------------------------Figures------------------------------------------

Figure.1<-plot()

#----------------------------------------References----------------------------------------
References<-stringr:str_c("
"\\subsection{KEGG}",
"\\bibitem[400]{key4000} Kanehisa, Furumichi, M., Tanabe, M., Sato, Y., and Morishima, K.", 
 "\\newblock KEGG: new perspectives on genomes, pathways, diseases and drugs.", 
 "\\newblock Nucleic Acids Res. 45, D353-D361 (2017).",
 "",
 "\\bibitem[401]{key4001} Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.", 
 "\\newblock KEGG as a reference resource for gene and protein annotation.", 
 "\\newblock Nucleic Acids Res. 44, D457-D462 (2016).",
 "",
 "\\bibitem[402]{key4002} Kanehisa, M. and Goto, S.", 
 "\\newblock KEGG: Kyoto Encyclopedia of Genes and Genomes.", 
 "\\newblock Nucleic Acids Res. 28, 27-30 (2000). ",
 "",
 "\\bibitem[403]{key4003}Petri, V., Jayaraman, P., Tutaj, M., Hayman, G. T., Smith, J. R., De Pons, J., … Jacob, H. J. (2014).", 
 "\\newblock The pathway ontology – updates and applications.", 
 "\\newblock Journal of Biomedical Semantics, 5, 7. http://doi.org/10.1186/2041-1480-5-7")

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
#---------------------------------------Function Library-----------------------------------
