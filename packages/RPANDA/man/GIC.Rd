\name{GIC.fit_pl.rpanda}
\alias{GIC.fit_pl.rpanda}

\title{
Generalized Information Criterion (GIC) to compare models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\description{
The GIC allows comparing models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\usage{

\method{GIC}{fit_pl.rpanda}(object, ...)
  
}
\arguments{
  \item{object}{
An object of class "fit_pl.rpanda". See ?fit_t_pl}
  \item{...}{
Options to be passed through.}
}


\value{
		a list with the following components
  
  \item{LogLikelihood}{the log-likelihood estimated for the model with estimated parameters}
  \item{GIC}{the GIC criterion}
  \item{bias}{the value of the bias term estimated to compute the GIC}

}

\details{
\code{GIC} allows comparing the fit of various models estimated by Penalized Likelihood (see ?fit_t_pl). It's a wrapper to the \code{gic_criterion} function.
}

	
\references{
Konishi S., Kitagawa G. 1996. Generalised information criteria in model selection. Biometrika. 83:875-890.

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68: 93-116.
}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{sim_t_env}}  
\code{\link{gic_criterion}}, 
\code{\link{fit_t_pl}}
\code{\link{mvgls}}
}

\examples{
\donttest{
require(mvMORPH)
set.seed(1)
n <- 32 # number of species
p <- 40 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- Posdef(p)      # a random symmetric matrix (covariance)
# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

fit1 <- fit_t_pl(Y, tree, model="BM", method="RidgeAlt")
fit2 <- fit_t_pl(Y, tree, model="OU", method="RidgeAlt")

GIC(fit1); GIC(fit2)
}
}