\name{getTipDistribution}
\alias{getTipDistribution}
\title{
Distribution of tip trait values.
}
\description{
Computes the mean and variance of the tip trait distribution under a specified model of trait evolution.
}
\usage{
getTipDistribution(object, params, v)
}

\arguments{
\item{object}{
an object of class 'PhenotypicModel'
}
\item{params}{
vector of parameters, given in the same order as in the 'model' object.
}
\item{v}{
boolean specifying the verbose mode. Default value : FALSE.
}
}

\value{
\item{mean}{Expectation vector of the tip trait distribution.}
\item{Sigma}{Variance-covariance matrix of the tip trait distribution.}
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

#Creating a BM model
modelBM <- createModel(tree, 'BM')

#Tip trait distribution under the model :
getTipDistribution(modelBM, c(0,0,0,1))
}
