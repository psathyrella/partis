\name{modelSelection}
\alias{modelSelection}
\title{
Phenotypic model selection from tip trait data.
}
\description{
For each model taken as input, fits the model and returns its AIC value in a recap table.
}
\usage{
modelSelection(object, data)
}
\details{
Warning : This function relies on the standard R optimizer "optim". It may not always converge well. Please double check the convergence by trying distinct parameter sets for the initialisation.
}
\arguments{
\item{object}{a vector of objects of class 'PhenotypicModel'.}
\item{data}{vector of tip trait data.}
}

\value{
A recap table presenting the AIC value of each model.
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

