\name{simulateTipData}
\alias{simulateTipData}
\title{
Tip trait simulation under a model of phenotypic evolution.
}
\description{
Simulates tip trait data under a specified model of phenotypic evolution, with three distinct behaviours specified with the 'method' argument.
}
\usage{
simulateTipData(object, params, method, v)
}

\arguments{
\item{object}{an object of class 'PhenotypicModel'.}
\item{params}{vector of parameters, given in the same order as in the 'model' object.}
\item{method}{an integer specifying the behaviour of the function. If method = 1 (default value), the tip distribution is first computed, before returning a simulated dataset drawn in this distribution. If method = 2, the whole trajectory is simulated step by step, plotted, and returned. Otherwise, the whole trajectory is simulated step by step, and then returned without being plotted.}
\item{v}{boolean specifying the verbose mode. Default value : FALSE.}
}

\value{
a vector of trait values at the tips of the tree.
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

#Simulating tip traits under both models with distinct behaviours of the functions :
dataBM <- simulateTipData(modelBM, c(0,0,0,1))
dataOU <- simulateTipData(modelOU, c(0,0,1,5,1), method=1)
dataBM2 <- simulateTipData(modelBM, c(0,0,0,1), method=2)
}
