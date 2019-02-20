\name{getDataLikelihood}
\alias{getDataLikelihood}
\title{
Likelihood of tip trait values.
}
\description{
Computes -log( likelihood ) of tip trait data under a given set of parameters, and for a specified model of trait evolution.
}
\usage{
getDataLikelihood(object, data, params, v)
}

\arguments{
\item{object}{an object of class 'PhenotypicModel'.}
\item{data}{vector of tip trait data.}
\item{params}{vector of parameters, given in the same order as in the 'model' object.}
\item{v}{boolean specifying the verbose mode. Default value : FALSE.}
}

\value{
A numerical value : -log( likelihood ) of the model.
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

#Simulating tip traits under the model :
dataBM <- simulateTipData(modelBM, c(0,0,0,1))

#Likelihood of the data :
getDataLikelihood(modelBM, dataBM, c(0,0,0,1))
}
