\name{sim.convergence.geo}
\alias{sim.convergence.geo}

\title{
Simulation of trait data under the model of convergent character displacement described in Drury et al. 2017
}

\description{
Simulates the evolution of a continuous character that evolves depending on pairwise similarity in another, OU-evolving trait (e.g., a trait that covaries with resource use). sig2 and z0 are shared between two traits, max and alpha are for focal trait, OU parameters for non-focal trait
}

\usage{

sim.convergence.geo(phylo,pars, Nsegments=2500, plot=FALSE, geo.object)

}

\arguments{
  \item{phylo}{
	an object of type 'phylo' (see ape documentation)
  }

  \item{pars}{ 
  A matrix with a number of rows corresponding to the desired number of simulations, columns containing values for \eqn{sig2} in \code{[,1]}, \eqn{m} in \code{[,2]}, \eqn{alpha} in \code{[,3]}, \code{root.value} in \code{[,4]}, \eqn{psi} of the OU model for the non-focal, resource use trait in \code{[,5]}, and \eqn{theta} in the OU model for the non-focal resource use trait in \code{[,6]} 
  }

\item{Nsegments}{
the minimum number of time steps to simulate
}

\item{plot}{
if \code{TRUE}, returns two plots: the top plot is focal trait undergoing convergence, the bottom plot is non-focal trait evolving under BM or OU
}

\item{geo.object}{
geography object created using \link{CreateGeoObject}
}
}



\value{

A list of two matrices with the simulated values for each lineage (one simulation per row; columns correspond to species) for trait1 (focul trait undergoing convergence) and non.focal (resource-use trait that determines strength of convergence in trait1)

}

\details{
Adjusting \code{Nsegments} will impact the length of time the simulations take. 
The length of each segment (\code{max(nodeHeights(phylo))/Nsegments}) 
should be much smaller than the smallest branch (\code{min(phylo$edge.length)}).
}

\references{
	Drury, J., Grether, G., Garland Jr., T., and Morlon, H. 2017. A review of phylogenetic methods for assessing the influence of interspecific interactions on phenotypic evolution. \emph{Systematic Biology} 
}

\author{
	J.P. Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}
}

\seealso{
\code{\link{CreateGeoObject}}
}

\examples{
\donttest{

data(Anolis.data)
phylo<-Anolis.data$phylo
geo.object<-Anolis.data$geography.object

#simulate with the OU process present and absent
pars<-expand.grid(0.05,-0.1,1,0,c(2,0),0)
sim.convergence.geo(phylo,pars,Nsegments=2500, plot=FALSE, geo.object)
}

}

