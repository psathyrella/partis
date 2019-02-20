\name{likelihood_t_DD_geog}
\alias{likelihood_t_DD_geog}

\title{
Likelihood of a dataset under diversity-dependent models with biogeography.
}
\description{
Computes the likelihood of a dataset under either the linear or exponential diversity dependent model with specified \code{sigma2} and slope values and with a \code{geography.object} formed using \code{\link{CreateGeoObject}}.  
}

\usage{
likelihood_t_DD_geog(phylo, data, par,geo.object,model=c("DDlin","DDexp"),maxN=NA)
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  
  \item{data}{
a named vector of continuous data with names corresponding to \code{phylo$tip.label}
}

 \item{par}{
a vector listing a value for \code{log(sig2)} (see Note) and either \code{b} (for the linear diversity dependent model) or \code{r} (for the exponential diversity dependent model), in that order.
}

  \item{geo.object}{
a list of sympatry through time created using \code{\link{CreateGeoObject}}
}

 \item{model}{
model chosen to fit trait data, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}

 \item{maxN}{
when fitting \code{DDlin} model, it is necessary to specify the maximum number of sympatric lineages to ensure that the rate returned does not correspond to negative sig2 values at any point in time (see Details).
}

}

\details{
	When specifying \code{par}, \code{log(sig2)} (see Note) must be listed before the slope parameter (\code{b} or \code{r}).

	maxN can be calculated using \code{maxN=max(vapply(geo.object$geography.object,function(x)max(rowSums(x)),1))}, where geo.object is the output of \code{\link{CreateGeoObject}}

}

\note{
	To stabilize optimization, this function exponentiates the input \code{sig2} value, thus the user must input the log(sig2) value to compute the correct log likelihood (see example).
}	

\value{
  the negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny,  \code{sig2} and slope values, and \code{geography.object}.
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

		Julien Clavel
}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{CreateGeoObject}}
\code{\link{likelihood_t_DD}}

}


\examples{
data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data
geography.object <- Anolis.data$geography.object

# Compute the likelihood with geography using ML parameters for fit without geography
par <- c(log(0.01153294),-0.0006692378)
maxN<-max(vapply(geography.object$geography.object,function(x)max(rowSums(x)),1))
lh <- -likelihood_t_DD_geog(phylo,pPC1,par,geography.object,model="DDlin",maxN=maxN)
}
