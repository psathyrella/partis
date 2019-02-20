\name{likelihood_sgd}
\alias{likelihood_sgd}
\title{
Likelihood of a phylogeny under the SGD model
}
\description{
Computes the likelihood of a phylogeny under the SGD model with exponential increasing of the metacommunity, and potentially missing extant species. Notations follow Manceau et al. (2015).
}
\usage{
likelihood_sgd(phylo, tot_time, b, d, nu, f)
}

\arguments{
\item{phylo}{
an object of type 'phylo' (see ape documentation)
}
\item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
\item{b}{
the (constant) birth rate of individuals in the model.
}
\item{d}{
the (constant) death rate of individuals in the model.
}
\item{nu}{
the (constant) mutation rate of individuals in the model.
}
\item{f}{
the fraction of extant species included in the phylogeny
}
}


\value{
the likelihood value of the phylogeny, given the model and the parameter values b, d, nu.
}

\references{
Manceau M., Lambert A., Morlon H. (2015) Phylogenies support out-of-equilibrium models of biodiversity Ecology Letters 18: 347-356
}

\author{
M Manceau
}

\examples{
data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
b <- 1e6
d <- 1e6-0.5
nu <- 0.6
f <- 87/89
#lh <- likelihood_sgd(Cetacea, tot_time, b, d, nu, f)
}
