\name{PhenotypicACDC-class}
\Rdversion{1.1}
\docType{class}
\alias{PhenotypicACDC-class}

\title{Class \code{"PhenotypicACDC"}}
\description{
Subclass of the PhenotypicModel class intended to represent the model of ACcelerating or DeCelerating phenotypic evolution.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PhenotypicACDC", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{matrixCoalescenceTimes}:}{Object of class \code{"matrix"} ~~ }
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
\section{Extends}{
Class \code{"\linkS4class{PhenotypicModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getTipDistribution}{\code{signature(object = "PhenotypicACDC")}: ... }
	 }
}
\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology, and the associated Supplementary material.
}
\author{
Marc Manceau
}

\examples{
showClass("PhenotypicACDC")
}
\keyword{classes}
