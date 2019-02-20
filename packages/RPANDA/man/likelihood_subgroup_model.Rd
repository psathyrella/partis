\name{likelihood_subgroup_model}
\alias{likelihood_subgroup_model}

\title{
Likelihood of a dataset under models with biogeography fit to a subgroup.
}
\description{
Computes the likelihood of a dataset under either the linear or exponential diversity dependent model with specified \code{sigma2} and slope values and with a \code{geography.object} formed using \code{\link{CreateGeoObject}}.  
}

\usage{
likelihood_subgroup_model(data,phylo,geography.object,model=c("MC","DDexp","DDlin"),
	par,return.z0=FALSE,maxN=NULL)

}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation) produced as "map" from CreateGeobyClassObject. NB: the length of this object need not match number of items in data, since map may include tips outside of group with some part of their branch in the group
}
  
  \item{data}{
a named vector of continuous data for a subgroup of interest with names corresponding to \code{phylo$tip.label}
}

  \item{geography.object}{
a list of sympatry/group membership through time created using \code{\link{CreateGeobyClassObject}}
}

 \item{model}{
model chosen to fit trait data, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}

 \item{par}{
a vector listing a value for \code{log(sig2)} (see Note) and either \code{b} (for the linear diversity dependent model) or \code{r} (for the exponential diversity dependent model), in that order.
}

 \item{return.z0}{
logical indicating whether to return an estimate of the trait value at the root given the parameter values (if \code{TRUE}, function returns root value rather than negative log-likelihood)
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
  The negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny,  \code{sig2} and slope values, and \code{geography.object}.
  
  If \code{return.z0=TRUE}, the estimated root value for the par values is returned instead of the negative log-likelihood.
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
data(BGB.examples)


Canidae.phylo<-BGB.examples$Canidae.phylo
dummy.group<-c(rep("B",3),rep("A",12),rep("B",2),rep("A",6),rep("B",5),rep("A",6))
names(dummy.group)<-Canidae.phylo$tip.label

\donttest{

Canidae.simmap<-make.simmap(Canidae.phylo, dummy.group)

set.seed(123)
Canidae.data<-rnorm(length(Canidae.phylo$tip.label))
names(Canidae.data)<-Canidae.phylo$tip.label
Canidae.A<-Canidae.data[which(dummy.group=="A")]
Canidae.geobyclass.object<-CreateGeobyClassObject(phylo=Canidae.phylo, 
	simmap=Canidae.simmap, trim.class="A", ana.events=BGB.examples$Canidae.ana.events, 
	clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE, rnd=5)

par <- c(log(0.01),-0.0005)
maxN<-max(vapply(Canidae.geobyclass.object$geography.object$geography.object, 
	function(x)max(rowSums(x)),1))

lh <- -likelihood_subgroup_model(data=Canidae.A, phylo=Canidae.geobyclass.object$map, 
	geo.object=Canidae.geobyclass.object$geography.object, model="DDlin", par=par, 
	return.z0=FALSE, maxN=maxN)
	
	}

}
