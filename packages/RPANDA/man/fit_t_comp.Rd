\name{fit_t_comp}
\alias{fit_t_comp}

\title{
Fits models of trait evolution incorporating competitive interactions
}
\description{
Fits matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution to a given dataset and phylogeny.
}
\usage{
fit_t_comp(phylo, data, model=c("MC","DDexp","DDlin"), pars=NULL, geography.object=NULL)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{data}{
a named vector of trait values with names matching \code{phylo$tip.label}
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model of Nuismer & Harmon 2014, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}
  \item{geography.object}{
a list of sympatry through time created using \code{\link{CreateGeoObject}}
}
  \item{pars}{
vector specifying starting parameter values for maximum likelihood optimization. If unspecified, default values are used (see Details)
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
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

		Julien Clavel 
}

\seealso{
\code{\link{sim_t_comp}}
\code{\link{CreateGeoObject}}
\code{\link{likelihood_t_MC}}
\code{\link{likelihood_t_MC_geog}}
\code{\link{likelihood_t_DD}}
\code{\link{likelihood_t_DD_geog}}
}

\note{
In current version, the \code{S} parameter is restricted to take on negative values in MC + geography ML optimization.
}
\examples{

data(Anolis.data)
geography.object<-Anolis.data$geography.object
pPC1<-Anolis.data$data
phylo<-Anolis.data$phylo

\donttest{
#Fit three models without biogeography to pPC1 data
MC.fit<-fit_t_comp(phylo, pPC1, model="MC")
DDlin.fit<-fit_t_comp(phylo, pPC1, model="DDlin")
DDexp.fit<-fit_t_comp(phylo, pPC1, model="DDexp")

#Now fit models that incorporate biogeography, NOTE these models take longer to fit
MC.geo.fit<-fit_t_comp(phylo, pPC1, model="MC", geography.object=geography.object)
DDlin.geo.fit<-fit_t_comp(phylo, pPC1,model="DDlin", geography.object=geography.object)
DDexp.geo.fit<-fit_t_comp(phylo, pPC1, model="DDexp", geography.object=geography.object)
}

}
