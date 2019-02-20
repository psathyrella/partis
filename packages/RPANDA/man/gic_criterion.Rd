\name{gic_criterion}
\alias{gic_criterion}

\title{
Generalized Information Criterion (GIC) to compare models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\description{
The GIC allows comparing models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\usage{


gic_criterion(Y, tree, model="BM", method=c("RidgeAlt", "RidgeArch", "LASSO", "ML", 
				"RidgeAltapprox", "LASSOapprox"), targM=c("null", 
				"Variance", "unitVariance"), param=NULL, 
				tuning=0, REML=TRUE, ...)
  
  
}
\arguments{
  \item{Y}{
A matrix of phenotypic traits values (the variables are represented as columns)
}
  \item{tree}{
An object of class 'phylo' (see ape documentation)
}
  \item{model}{
The evolutionary model, "BM" is Brownian Motion, "OU" is Ornstein-Uhlenbeck, "EB" is Early Burst, and "lambda" is Pagel's lambda transformation.
}
  \item{method}{
The penalty method. "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator, "ML": Maximum Likelihood.
}
  \item{targM}{
The target matrix used for the Ridge regularizations. "null" is a null target, "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt" methods.
}

 \item{param}{
Parameter for the evolutionary model (see "model" above).
}

 \item{tuning}{
The tuning/regularization parameter.
}

  \item{REML}{
Use REML (default) or ML for estimating the parameters.
}
\item{...}{
Additional options. Not used yet.}
}


\value{
		a list with the following components
  
  \item{LogLikelihood}{the log-likelihood estimated for the model with estimated parameters}
  \item{GIC}{the GIC criterion}
  \item{bias}{the value of the bias term estimated to compute the GIC}

}

\details{
\code{gic_criterion} allows comparing the fit of various models estimated by Penalized Likelihood (see ?fit_t_pl). Use the wrapper \code{GIC} instead for models fit with \code{fit_t_pl}.
}

\note{
	The tuning parameter is assumed to be zero when using the "ML" method.}
	
\references{
Konishi S., Kitagawa G. 1996. Generalised information criteria in model selection. Biometrika. 83:875-890.

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68: 93-116.
}

\author{J. Clavel}

\seealso{
\code{\link{GIC.fit_pl.rpanda}}, 
\code{\link{fit_t_pl}}
}

\examples{
\donttest{
require(mvMORPH)
set.seed(123)
n <- 32 # number of species
p <- 2 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- Posdef(p)      # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

# Compute the GIC for ML
gic_criterion(Y, tree, model="BM", method="ML", tuning=0) # ML

# Compare with PL?
#test <- fit_t_pl(Y, tree, model="BM", method="RidgeAlt")
#GIC(test)
}
}
