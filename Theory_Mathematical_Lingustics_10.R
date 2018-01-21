library(ontologyIndex)
library(ontologyPlot)
library(ontologySimilarity)
library(stringr)
library(stringi)
library(readr)
library(readxl)
library(xtable)

#------------------------------------------Data---------------------------------------------

Category.A<-list(A=character(0),A.1="A",A.2="A.1",A.3="A")
Category.B<-list(B=character(0), B.1="B",B.2="B.1",B.3="B.1")

#------------------------------------------Ontology Enrichment---------------------------------------

A.ontology <- ontology_index(parents=Category.A)
B.ontology <- ontology_index(parents=Category.B)
unclass(A.ontology)
unclass(B.ontology)

A.ontology.frequencies<-get_term_frequencies(A.ontology,list('A','A.1'))
A.ontology.info<-get_term_info_content(A.ontology,list('A'))
B.ontology.frequencies<-get_term_frequencies(B.ontology,list('B.1'))
B.ontology.terms <- remove_links(B.ontology, get_ancestors(B.ontology, c("B.1")))

print(A.ontology)
A.ontology.descendents<-get_descendants(A.ontology,roots='A')
A.ontology.descendents.intersection<-intersection_with_descendants(A.ontology,c('A','B'),c('A'))

#-----------------------------------------Group Similarity------------------------------------

A.ontology.information.content<- descendants_IC(A.ontology)
A.ontology.term.sets <- replicate(simplify=TRUE, n=8, expr=minimal_set(A.ontology, sample(A.ontology$id, size=3)))
A.ontology.sim.mat <- get_sim_grid(ontology=A.ontology, term_sets=A.ontology.term.sets)
A.ontology.dist.mat <- max(A.ontology.sim.mat) - A.ontology.sim.mat

group <- 1:3
A.ontology.similarity<-get_sim_p_from_ontology(ontology=A.ontology, information_content=A.ontology.information.content,
                        term_sets=A.ontology.term.sets,group=group)
A.ontology.similarity.pvalue<-get_sim_p(A.ontology.sim.mat,group=group)

A.ontology.group.sim <- get_sim(A.ontology.sim.mat, group=group)
A.ontology.sim.mat.samples <- sample_group_sim(A.ontology.sim.mat, group_size=length(group))

Ontology.information.content.df<-data.frame()
Ontology.information.content.df<-rbind(c(descendants_IC(A.ontology)),
                                         c(descendants_IC(B.ontology)))
rownames(Ontology.information.content.df)<-c("A","B")
colnames(Ontology.information.content.df)<-c("Ontology Root","1","2","3")
#------------------------------------------Tables-------------------------------------------

Table.1<-xtable(Ontology.information.content.df)

#------------------------------------------Figures-----------------------------------------
Figure.1<-onto_plot(A.ontology, terms=A.ontology.descendents)
frequencies <-B.ontology.frequencies
names(frequencies) <- B.ontology.terms
Figure.2<-onto_plot(B.ontology, term_sets=B.ontology.terms, frequencies=frequencies, 
                    fillcolor=colour_by_frequency)
Figure.3<-plot(hclust(as.dist(A.ontology.dist.mat)))
Figure.4<-hist(A.ontology.sim.mat.samples)
abline(v=A.ontology.group.sim, col="blue")
#-----------------------------------------References---------------------------------------

Reference.1<-c("Smith, Cynthia L., Carroll-Ann W. Goldsmith, and Janan T. Eppig.", 
               "The Mammalian Phenotype Ontology as a tool for annotating, analyzing and comparing phenotypic information.",
               "Genome biology 6.1 (2004): 1.")

#----------------------------------------Function Library----------------------------------
