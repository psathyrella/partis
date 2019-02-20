\name{fit_sgd}
\alias{fit_sgd}
\title{
Maximum likelihood fit of the SGD model
}

\description{
Fits the SGD model with exponential growth of the metacommunity, by maximum likelihood. Notations follow Manceau et al. (2015)
}

\usage{
fit_sgd(phylo, tot_time, par, f=1, meth = "Nelder-Mead")
}

\arguments{
\item{phylo}{
an object of type 'phylo' (see ape documentation)
}
\item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages)
}
\item{par}{
a numeric vector of initial values for the parameters (b,d,nu) to be estimated (these values are used by the optimization algorithm)
}
\item{f}{
the fraction of extant species included in the phylogeny
}
\item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
}

\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{par}{a numeric vector of estimated values of b (birth), b-d (growth) and nu (mutation)}
}

\note{While b-d and nu can in general be well estimated, the likelihood surface is quite flat whith respect to b, such that the estimated b can vary a lot depending on the choice of the initial parameter values. Estimates of b should not be trusted.
	}

\author{
M Manceau
}

\references{
Manceau, M., Lambert, A., Morlon, H. (2015) Phylogenies support out-of-equilibrium models of biodiversity Ecology Letters 18: 347-356
}

\seealso{
\code{\link{likelihood_sgd}}
}

\examples{
# Some examples may take a little bit of time. Be patient!
data(Calomys)
tot_time <- max(node.age(Calomys)$ages)
par_init <- c(1e7, 1e7-0.5, 1)
#fit_sgd(Calomys, tot_time, par_init, f=11/13)
}
