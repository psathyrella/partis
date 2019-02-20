\name{CreateClassObject}
\alias{CreateClassObject}

\title{
Create class object
}
\description{
This function returns names of internode intervals, named descendants of each node,
and a class object formatted in a way that can be
passed to \code{\link{CreateGeobyClassObject}}

}
\usage{

CreateClassObject(simmap,rnd=5)

}

\arguments{

  \item{simmap}{
  stochastic map from \code{\link{make.simmap}} in \code{phytools}
}
  \item{rnd}{ 
  integer indicating the number of decimal places to which times should be rounded (default value is 5) (see \code{\link{round}})
}

}

\value{
	a list with the following components:
  \item{class.object}{a list of matrices specifying the state of each branch during each internode interval (see Details)}
  \item{times}{a vector containing the time since the root of the tree at which nodes or changes in biogeography occur (used internally in other functions)}
  \item{spans}{a vector specifying the distances between times (used internally in other functions)}
}

\details{
This function formats the class object so that it can be correctly
passed to the numerical integration performed in \code{\link{fit_t_comp_subgroup}}. 

}

\references{

	Drury, J., Tobias, J., Burns, K., Mason, N., Shultz, A., and Morlon, H. in review. Contrasting impacts of competition on ecological and social trait evolution in songbirds. \emph{PLOS Biology}.

	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}}

\seealso{
\code{\link{fit_t_comp_subgroup}},\code{\link{CreateGeobyClassObject}} 
}

\examples{

data(Anolis.data)

#Create a make.simmap object
require(phytools)
geo<-c(rep("cuba",7),rep("hispaniola",9),"puerto_rico")
names(geo)<-Anolis.data$phylo$tip.label
stochastic.map<-phytools::make.simmap(Anolis.data$phylo, 
									geo, model="ER", nsim=1)
CreateClassObject(stochastic.map)

}
