\name{likelihood_t_MC}
\alias{likelihood_t_MC}

\title{
Likelihood of a dataset under the matching competition model.
}
\description{
Computes the likelihood of a dataset under the matching competition model with specified \code{sigma2} and \code{S} values. 
}
\usage{
likelihood_t_MC(phylo, data, par)
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  
  \item{data}{
a named vector of continuous data with names corresponding to \code{phylo$tip.label}
}
 \item{par}{
a vector listing a value for \code{log(sig2)} (see Note) and \code{S} (parameters of the matching competition model), in that order
}
}

\details{
	When specifying \code{par}, \code{log(sig2)} must be listed before \code{S}.
}

\note{
	To stabilize optimization, this function exponentiates the input \code{sig2} value, thus the user must input the log(sig2) value to compute the correct log likelihood (see example).
}	

\value{
  the negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny and \code{sig2} and \code{S} values
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

		Julien Clavel
}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{likelihood_t_MC_geog}}

}


\examples{
data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data

# Compute the likelihood that the S value is twice the ML estimate
par <- c(0.0003139751, (2*-0.06387258))
lh <- -likelihood_t_MC(phylo,pPC1,par)
}
