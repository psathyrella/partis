plot_BICompare <- function(phylo,BICompare)
{
  if (!inherits(BICompare, "BICompare"))
      stop("object \"BICompare\" is not of class \"BICompare\"")

t<-max(BICompare[[2]])
col_edge<-rainbow(t)[BICompare[[2]][phylo$edge[,2]]]
col_tip<-rainbow(t)[BICompare[[2]][1:length(phylo$tip.label)]]
plot(phylo,edge.color=col_edge,tip.color=col_tip,type="fan",cex=0.4)

  
  }

