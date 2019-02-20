\name{sim.divergence.geo}
\alias{sim.divergence.geo}

\title{
Simulation of trait data under the model of divergent character displacement described in Drury et al. 2017
}
\description{
Simulates the evolution of a continuous character under a model of evolution where trait values are repelled according to between-species similarity in trait values, taking into account biogeography using a biogeo.object formatted from RPANDA (see CreateGeoObject function in RPANDA package)
}

\usage{

sim.divergence.geo(phylo,pars, Nsegments=2500, plot=FALSE, geo.object)

}
\arguments{

\item{phylo}{
	a phylogenetic tree
}
\item{pars}{
	A matrix with a number of rows corresponding to the desired number of simulations, columns containing values for \eqn{sig2} in \code{[,1]}, \eqn{m} in \code{[,2]}, \eqn{alpha} in \code{[,3]}, \code{root.value} in \code{[,4]}, \eqn{psi} of the OU model in \code{[,5]}, and \eqn{theta} in the OU model in \code{[,6]}
}
\item{Nsegments}{
the minimum number of time steps to simulate
}
\item{plot}{
logical indicating whether to plot the simulated trait values at each time step
}
\item{geo.object}{
geography object created using \link{CreateGeoObject}
}
}


\value{

A matrix with the simulated values for each lineage (one simulation per row; columns correspond to species) 

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
	F. Hartig
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
pars<-expand.grid(0.05,2,1,0,c(2,0),0)
sim.divergence.geo(phylo,pars,Nsegments=2500, plot=FALSE, geo.object)
}

}

