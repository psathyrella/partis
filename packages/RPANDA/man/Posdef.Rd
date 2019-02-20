\name{Posdef}
\alias{Posdef}

\title{
Positive definite symmetric matrices
}
\description{
Generates a positive definite and symmetric matrix with specified eigen-values
}
\usage{


Posdef(p, ev = rexp(p, 1/100))
  
  
}
\arguments{
  \item{p}{
The dimension of the matrix}

  \item{ev}{
The eigenvalues. If not specified, eigenvalues are taken from an exponential distribution.}

}


\value{
		Returns a symmetric positive-definite matrix with eigen-values = ev.

}

\details{
\code{Posdef} generates random positive definite covariance matrices with specified eigen-values that can be used to simulate multivariate datasets (see Uyeda et al. 2015 - and supplied R codes). 
}

	
\references{
Uyeda J.C., Caetano D.S., Pennell M.W. 2015. Comparative Analysis of Principal Components Can be Misleading. Syst. Biol. 64:677-689.
	

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68:93-116.
}

\author{J. Clavel}

\seealso{
\code{\link{GIC.fit_pl.rpanda}}, 
\code{\link{fit_t_pl}}
\code{\link{phyl.pca_pl}}
}

\examples{
\donttest{
require(mvMORPH)
set.seed(123)
n <- 32 # number of species
p <- 40 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- Posdef(p) # a random symmetric matrix (covariance)
# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

test <- fit_t_pl(Y, tree, model="BM", method="RidgeAlt")
GIC(test)

}
}