\name{sim_t_comp}
\alias{sim_t_comp}

\title{
Recursive simulation (root-to-tip) of competition models 
}
\description{
Simulates datasets for a given phylogeny under matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution. Simulations are carried out from the root to the tip of the tree.
}
\usage{

sim_t_comp(phylo,pars,root.value,Nsegments=1000,model="MC,DDexp,DDlin")

}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{pars}{
a vector containing the two parameters for the chosen model; all models require \code{sig2}, and additionally, the MC model requires \code{S},
 specifying the level of competition (larger negative values correspond to higher levels of competition), the DDlin model 
 requires \code{b} and DDexp require \code{r}, the slope parameters (negative in cases of decline in evolutionary 
 rates with increasing diversity). \code{sig2} must be listed first.
}
  \item{root.value}{
a number specifying the trait value for the ancestor
}
  \item{Nsegments}{
a value specifying the total number of time segments to simulate across for the phylogeny (see Details)
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model of Nuismer & Harmon 2014, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}
}

\value{
	a named vector with simulated trait values for \eqn{n} species in the phylogeny
}

\details{
Adjusting \code{Nsegments} will impact the length of time the simulations take. 
The length of each segment (\code{max(nodeHeights(phylo))/Nsegments}) 
should be much smaller than the smallest branch (\code{min(phylo$edge.length)}).
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{J Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}}

\seealso{
\code{\link{fit_t_comp}}
}

\examples{

data(Cetacea)

\donttest{
# Simulate data under the matching competition model
MC.data<-sim_t_comp(Cetacea,pars=c(sig2=0.01,S=-0.1),root.value=0,Nsegments=1000,model="MC")

# Simulate data under the diversity dependent linear model
DDlin.data<-sim_t_comp(Cetacea,pars=c(sig2=0.01,b=-0.0001),root.value=0,Nsegments=1000,
	model="DDlin")

# Simulate data under the diversity dependent linear model
DDexp.data<-sim_t_comp(Cetacea,pars=c(sig2=0.01,r=-0.01),root.value=0,Nsegments=1000,model="DDexp")
}

}
