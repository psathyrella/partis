\name{BGB.examples}
\alias{BGB.examples}

\title{
BioGeoBEARS stochastic maps
}
\description{
Phylogenies and example stochastic maps for Canidae (from an unstratified BioGeoBEARS analysis) and Ochotonidae (from a stratified BioGeoBEARS analysis)
}
\usage{
data(BGB.examples)
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Matzke, N. 2014. Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. \emph{Systematic Biology} 63:951-970.
	}
	
\seealso{\code{\link{CreateGeoObject_BioGeoBEARS}}}


\examples{
data(BGB.examples)
par(mfrow=c(1,2))
plot(BGB.examples$Canidae.phylo)
plot(BGB.examples$Ochotonidae.phylo)

}
