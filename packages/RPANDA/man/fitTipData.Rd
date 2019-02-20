\name{fitTipData}
\alias{fitTipData}
\title{
Maximum likelihood estimators of a model's parameters
}
\description{
Finds the maximum likelihood estimators of the parameters, returns the likelihood and the inferred parameters.
}
\usage{
fitTipData(object, data, params0, GLSstyle, v)
}
\details{
Warning : This function uses the standard R optimizer "optim". It may not always converge well. Please double check the convergence by trying distinct parameter sets for the initialisation.
}
\arguments{
\item{object}{an object of class 'PhenotypicModel'.}
\item{data}{vector of tip trait data.}
\item{params0}{vector of parameters used to initialize the optimization algorithm. Default value is NULL, in which case the optimization procedure starts with the vector 'params0' specified within the 'model' object. }
\item{GLSstyle}{boolean specifying the way the mean trait value at the root is estimated. Default value is FALSE in which case the mean at the root is considered as any other parameter. If TRUE, the mean value at the root is estimated with the GLS method, as explained, e.g. in Hansen 1997.}
\item{v}{boolean specifying the verbose mode. Default value : FALSE.}
}

\value{
\item{value}{A numerical value : the lowest -log( likelihood ) value found during the optimization procedure.}
\item{inferredParams}{The maximum likelihood estimators of the model's parameters.}
\item{convergence}{An integer code specifying the convergence of the optim function. Please refer to the optim function help files.}
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

#Fitting the model to the data
fitTipData(modelBM, dataBM, v=TRUE)
}
