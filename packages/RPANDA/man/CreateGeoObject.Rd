\name{CreateGeoObject}
\alias{CreateGeoObject}

\title{
Create biogeography object
}
\description{
This function returns names of internode intervals, named descendants of each node,
and a geography object formatted in a way that can be
passed to \code{\link{fit_t_comp}}

}
\usage{

CreateGeoObject(phylo,map)

}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{map}{
either a matrix modified from \code{phylo$edge} or a phylo object created using \code{\link{make.simmap}} (see Details and Examples)
}
}
\value{
	a list with the following components:
  \item{geography.object}{a list of matrices specifying sympatry (1) or allopatry (0) for each species pair for each internode interval (see Details)}
  \item{times}{a vector containing the time since the root of the tree at which nodes or changes in biogeography occur (used internally in other functions)}
  \item{spans}{a vector specifying the distances between times (used internally in other functions)}
}

\details{
This function should be used to format the geography object so that it can be correctly
passed to the numerical integration performed in \code{\link{fit_t_comp}}. 

The \code{map} can either be a matrix formed by specifying the region in which each branch specified by \code{phylo$edge} existed, or a stochastic map stored as a phylo object output from \code{\link{make.simmap}} (see Examples).
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}}

\seealso{
\code{\link{fit_t_comp}} 
}

\examples{

data(Anolis.data)
#Create a geography.object with a modified edge matrix
#First, specify which region each branch belonged to:
Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
Anolis.map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
CreateGeoObject(Anolis.data$phylo,map=Anolis.map)

#Create a geography.object with a make.simmap object
#First, specify which region each branch belonged to:
require(phytools)
geo<-c(rep("cuba",7),rep("hispaniola",9),"puerto_rico")
names(geo)<-Anolis.data$phylo$tip.label
stochastic.map<-phytools::make.simmap(Anolis.data$phylo, 
							geo, model="ER", nsim=1)
CreateGeoObject(Anolis.data$phylo,map=stochastic.map)

}
