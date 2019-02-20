\name{PhenotypicModel-class}
\Rdversion{1.1}
\docType{class}
\alias{PhenotypicModel-class}
\alias{[<-,PhenotypicModel,ANY,ANY,ANY-method}
\alias{[,PhenotypicModel,ANY,ANY,ANY-method}
\alias{print,PhenotypicModel-method}
\alias{show,PhenotypicModel-method}

\title{Class \code{"PhenotypicModel"}}
\description{
This class describes a model of phenotypic evolution running on a phylogenetic tree, with or without interactions between lineages.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PhenotypicModel", ...)}.
Alternatively, you may just want to use the "createModel" function for predefined models.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{period}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{aAGamma}:}{Object of class \code{"function"} ~~ }
    \item{\code{numbersCopy}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{numbersPaste}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{initialCondition}:}{Object of class \code{"function"} ~~ }
    \item{\code{paramsNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{constraints}:}{Object of class \code{"function"} ~~ }
    \item{\code{params0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{tipLabels}:}{Object of class \code{"character"} ~~ }
    \item{\code{tipLabelsSimu}:}{Object of class \code{"character"} ~~ }
    \item{\code{comment}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "PhenotypicModel", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{[}{\code{signature(x = "PhenotypicModel", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{fitTipData}{\code{signature(object = "PhenotypicModel")}: ... }
    \item{getDataLikelihood}{\code{signature(object = "PhenotypicModel")}: ... }
    \item{getTipDistribution}{\code{signature(object = "PhenotypicModel")}: ... }
    \item{modelSelection}{\code{signature(object = "PhenotypicModel")}: ... }
    \item{print}{\code{signature(x = "PhenotypicModel")}: ... }
    \item{show}{\code{signature(object = "PhenotypicModel")}: ... }
    \item{simulateTipData}{\code{signature(object = "PhenotypicModel")}: ... }
	 }
}
\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology, and the associated Supplementary material.
}
\author{
Marc Manceau
}

\examples{
showClass("PhenotypicModel")
}
\keyword{classes}
