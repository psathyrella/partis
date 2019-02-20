\name{sim_env_bd}
\alias{sim_env_bd}

\title{
Simulate birth-death tree dependent on an environmental curve
}
\description{
Simulates a birth-death tree (starting with one lineage) with speciation and/or extinction rate that varies as a function of an input environmental curve. Notations follow Morlon et al. PNAS 2011 and Condamine et al. ELE 2013.
}
\usage{
sim_env_bd(env_data, f.lamb, f.mu, lamb_par, mu_par, df=NULL, time.stop=0, 
			return.all.extinct=TRUE, prune.extinct=TRUE)
}

\arguments{

  \item{env_data}{
environmental data, given as a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{time.stop}{
the age of the phylogeny.
}
  \item{f.lamb}{
a function specifying the hypothesized functional form of the variation of the speciation rate \eqn{\lambda} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the third arguement is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
}
  \item{f.mu}{
a function specifying the hypothesized functional form of the variation of the extinction rate \eqn{\mu} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the second argument is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
}
  \item{lamb_par}{
a numeric vector of initial values for the parameters of f.lamb to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong.
}
  \item{mu_par}{
a numeric vector of initial values for the parameters of f.mu to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model without extinction (for example), mu_par should be empty (vector of length 0). Otherwise aic values will be wrong.
}
\item{df}{
 the degree of freedom to use to define the spline. As a default, smooth.spline(env_data[,1], env_data[,2])$df is used. See \emph{sm.spline} for details.  
}
\item{return.all.extinct}{
 return all extinction lineages in simulated tree.  
}
\item{prune.extinct}{
 prune extinct lineages in simulated tree.  
}
  }


\value{
	a list with the following components
  \item{tree}{the simulated tree with number tips}
  \item{times}{the times of speciation events starting from the past}
  \item{nblineages}{the labels of surviving lineages and total number of surviving lineages}
}

\details{In the f.lamb and f.mu functions, time runs from the present to the past.
}
\note{
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}
\references{
Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	
Condamine, F.L., Rolland, J., and Morlon, H. (2013) Macroevolutionary perspectives to environmental change, Eco Lett 16: 72-85
}
\author{E Lewitus and H Morlon}
\seealso{
\code{\link{fit_env}}, \code{\link{fit_bd}}
}
\examples{
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df
# Simulates a tree with lambda varying as an exponential function of temperature
# and mu fixed to 0 (no extinction).  Here t stands for time and x for temperature.
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
#result_exp <- sim_env_bd(InfTemp,f.lamb,f.mu,lamb_par,mu_par,time.stop=10)
}
