\name{fit_t_comp_subgroup}
\alias{fit_t_comp_subgroup}

\title{
Fits models of trait evolution incorporating competitive interactions, restricting competition to occur only between members of a subgroup
}
\description{
Fits matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution to a given dataset, phylogeny, and stochastic maps of both subgroup membership and biogeography.
}
\usage{

fit_t_comp_subgroup(full.phylo, ana.events, clado.events, stratified=FALSE, map,data, 
	trim.class, model=c("MC","DDexp","DDlin"), par=NULL, method="Nelder-Mead", 
	bounds=NULL)
	
}

\arguments{
  \item{full.phylo}{
an object of type 'phylo' (see ape documentation) containing all of the tips used to estimate ancestral biogeography in BioGeoBEARS
}

  \item{ana.events}{
the "ana.events" table produced in BioGeoBEARS that lists anagenetic events in the stochastic map
}
  \item{clado.events}{
the "clado.events" table produced in BioGeoBEARS that lists cladogenetic events in the stochastic map
}
  \item{stratified}{
logical indicating whether the stochastic map was built from a stratified analysis in BioGeoBEARS
}
  \item{map}{
a phylo object created using \code{\link{make.simmap}} in phytools that contains reconstructed subgroup membership
}
  \item{data}{
a named vector of trait values for subgroup members with names matching \code{full.phylo$tip.label}
}
  \item{trim.class}{
  subgroup whose members are competing
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model of Nuismer & Harmon 2014, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}
  \item{par}{
vector specifying starting parameter values for maximum likelihood optimization. If unspecified, default values are used (see Details)
}
  \item{method}{
optimization algorithm to use (see \code{\link{optim}})
}
  \item{bounds}{
(optional) list of bounds to pass to optimization algorithm (see details at \code{\link{optim}})
}

}

\value{
	a list with the following elements:
  \item{LH}{maximum log-likelihood value}
  \item{aic}{Akaike Information Criterion value}
  \item{aicc}{AIC value corrected for small sample size}
  \item{free.parameters}{number of free parameters from the model}
  \item{sig2}{maximum-likelihood estimate of \code{sig2} parameter}
  \item{S}{maximum-likelihood estimate of \code{S} parameter of matching competition model (see Note)}
  \item{b}{maximum-likelihood estimate of \code{b} parameter of linear diversity dependence model}
  \item{r}{maximum-likelihood estimate of \code{r} parameter of exponential diversity dependence model}
  \item{z0}{maximum-likelihood estimate of \code{z0}, the value at the root of the tree}
  \item{convergence}{convergence diagnostics from \code{\link{optim}} function (see optim documentation)}
}

\details{
If unspecified, \code{par} takes the default values of \code{var(data)/max(nodeHeights(phylo))} for sig2 and 0 for either \code{S} for the matching competition model, 
\code{b} for the linear diversity dependence model, or \code{r} for the exponential diversity dependence model. Values can be manually entered as a vector with the first element
equal to the desired starting value for sig2 and the second value equal to the desired starting value for either \code{S}, \code{b}, or \code{r}. Note: since likelihood optimization
uses sig rather than sig2, and since the starting value for is exponentiated to stabilize the likelihood search, if you input a \code{par} value, the first value specifying sig2 
should be the log(sqrt()) of the desired sig2 starting value.

}

\references{
	
	Drury, J., Tobias, J., Burns, K., Mason, N., Shultz, A., and Morlon, H. in review. Contrasting impacts of competition on ecological and social trait evolution in songbirds. \emph{PLOS Biology}.

	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

}

\seealso{
\code{\link{likelihood_subgroup_model}}
\code{\link{CreateGeobyClassObject}}

}

\note{
In current version, the \code{S} parameter is restricted to take on negative values in MC + geography ML optimization.
}
\examples{


\donttest{

data(BGB.examples)


Canidae.phylo<-BGB.examples$Canidae.phylo
dummy.group<-c(rep("B",3),rep("A",12),rep("B",2),rep("A",6),rep("B",5),rep("A",6))
names(dummy.group)<-Canidae.phylo$tip.label


Canidae.simmap<-phytools::make.simmap(Canidae.phylo,dummy.group)

set.seed(123)
Canidae.data<-rnorm(length(Canidae.phylo$tip.label))
names(Canidae.data)<-Canidae.phylo$tip.label
Canidae.A<-Canidae.data[which(dummy.group=="A")]

fitA<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events, 
clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE,map=Canidae.simmap, 
data=Canidae.A,trim.class="A",model="DDexp")
	
}

}
