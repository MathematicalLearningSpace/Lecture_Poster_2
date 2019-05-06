%-----------------------------------------------------------------------%
%---------------------Classroom R Program Example-----------------------%
%-----------------------------------------------------------------------%




Module.Cell.Cycle.Model<-function(X,Y,Visualization=TRUE)
{
  library(KEGG.db);library(KEGGgraph);library(KEGGprofile);library(KEGGREST);library(PearsonDS)
  data(pro_pho_expr);data(pho_sites_count)
  Y<-"hsa04151.xml"
  path.id<-"04151"
  Y.30<-xml2::read_xml(Y)
  row.names(X1)
  groups <- rep("Expression", ncol(X1))
  groups <- factor(groups, levels = unique(groups), ordered = F)
  XML2database<-parse_XMLfile(pathway_id=path.id,species="hsa")
  Genes.df<-as.data.frame(XML2database)
  Gene.Expression <- intersect(Genes.df[, 1], row.names(X))
  X1<-pro_pho_expr
  Z<-X1[Gene.Expression,1:6]
 if(Visualization)
  {
    op <- par(mfrow = c(2,3),mar=c(3,3,3,3))
    hist(pro_pho_expr$Proteome.G1.phase,main="G1.Phase")
    hist(pro_pho_expr$Proteome.G1.S,main="G1.S")
    hist(pro_pho_expr$Proteome.Early.S,main="Early.S")
    hist(pro_pho_expr$Proteome.Late.S,main="Late.S")
    hist(pro_pho_expr$Proteome.G2.phase,main="G2.Phase")
    hist(pro_pho_expr$Proteome.Mitosis,main="Mitosis")
    par(op)
    op <- par(mfrow = c(2,3),mar=c(3,3,3,3))
    hist(Z[,1],main="G1.Phase")
    hist(Z[,2],main="G1.S")
    hist(Z[,3],main="Early.S")
    hist(Z[,4],main="Late.S")
    hist(Z[,5],main="G2.Phase")
    hist(Z[,6],main="Mitosis")
    par(op)
  }
  
  output<-list()
  return(output)
}
test.Module.Cell.Cycle.Model<-Module.Cell.Cycle.Model("1",TRUE)
test.Module.Cell.Cycle.Model
