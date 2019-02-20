\name{sim_sgd}
\alias{sim_sgd}
\title{
Algorithm for simulating a phylogenetic tree under the SGD model
}
\description{
Simulates a phylogeny arising from the SGD model with exponentially increasing metapopulation size. Notations follow Manceau et al. (2015).
}
\usage{
sim_sgd(tau, b, d, nu)
}

\arguments{
\item{tau}{
the simulation time, which corresponds to the length of the phylogeny
}
\item{b}{
the (constant) per-individual birth rate
}
\item{d}{
the (constant) per-individual death rate
}
\item{nu}{
the (constant) per-individual mutation rate
}
}

\value{
a phylogenetic tree of class "phylo" (see ape documentation)
}

\references{
Manceau M., Lambert A., Morlon H. (2015) Phylogenies support out-of-equilibrium models of biodiversity Ecology Letters 18: 347-356
}

\author{
M Manceau
}

\examples{
tau <- 10
b <- 1e6
d <- b-0.5
nu <- 0.6
tree <- sim_sgd(tau,b,d,nu)
plot(tree)
}
