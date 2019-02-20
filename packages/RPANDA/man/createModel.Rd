\name{createModel}
\alias{createModel}
\title{
Creation of a PhenotypicModel
}
\description{
Creates an object of class PhenotypicModel, intended to represent a model of trait evolution on a specific tree. DIstinct keywords correspond to different models, using one phylogenetic tree.
}
\usage{
createModel(tree, keyword)
}

\arguments{
\item{tree}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{keyword}{
a string specifying the model. Available models include "BM", "BM_from0", "BM_from0_driftless", "OU", "OU_from0", "ACDC", "DD", "PM", "PM_OUless".
}
}

\value{
the object of class "PhenotypicModel".
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

\examples{
#Loading an example tree
newick <- "((((A:1,B:0.5):2,(C:3,D:2.5):1):6,E:10.25):2,(F:6.5,G:8.25):3):1;"
tree <- read.tree(text=newick)

#Creating the models
modelBM <- createModel(tree, 'BM')
modelOU <- createModel(tree, 'OU')

#Printing basic or full informations on the model definitions
show(modelBM)
print(modelOU)
}
