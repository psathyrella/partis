\name{ancestral}
\alias{ancestral}

\title{
Estimation of traits ancestral states.
}
\description{
Reconstruct the ancestral states at the root (and possibly for each nodes) of a phylogenetic tree from models fit obtained using the \code{fit_t_XX} functions.}
\usage{

ancestral(object)
  
}
\arguments{
  \item{object}{
A model fit object obtained by the \code{fit_t_XX} class of functions.
}
}
 

\value{
	a list with the following components
  
  \item{root}{the reconstructed ancestral states at the root}
  \item{nodes}{the reconstructed ancestral states at each nodes (not yet implemented for all the methods)}

  }

\details{
\code{ancestral} is an S3 method that reconstruct the ancestral states at the root and possibly for each nodes of a phylogenetic tree from the models fit obtained by the \code{fit_t_XX} class of functions (e.g., \code{fit_t_pl}, \code{fit_t_comp} and \code{fit_t_env}). Ancestral states are estimated using generalized least squares (GLS; Martins & Hansen 1997, Cunningham et al. 1998 ).}

\note{
	The function is used internally in \code{phyl.pca_pl}.}
\references{

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68: 93-116.

Cunningham C.W., Omland K.E., Oakley T.H. 1998. Reconstructing ancestral character states: a critical reappraisal. Trends Ecol. Evol. 13:361-366.

Martins E.P., Hansen T.F. 1997. Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data. Am. Nat. 149:646-667.
}

\author{J. Clavel}

\seealso{
\code{\link{fit_t_pl}},
\code{\link{fit_t_env}},
\code{\link{phyl.pca_pl}},
\code{\link{GIC}}, 
\code{\link{gic_criterion}}
}

\examples{
\donttest{
require(mvMORPH)
set.seed(1)
n <- 32 # number of species
p <- 31 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- Posdef(p)      # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

# fit a multivariate BM with Penalized likelihood
fit <- fit_t_pl(Y, tree, model="BM", method="RidgeAlt")

# Perform the ancestral states reconstruction
anc <- ancestral(fit)

# retrieve the scores
head(anc$nodes)
}
}
