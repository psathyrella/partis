\name{phyl.pca_pl}
\alias{phyl.pca_pl}

\title{
Regularized Phylogenetic Principal Component Analysis (PCA).
}
\description{
Performs a principal component analysis (PCA) on a regularized evolutionary variance-covariance matrix obtained using the \code{fit_t_pl} function.}
\usage{

phyl.pca_pl(object, plot=TRUE, ...)
  
}
\arguments{
  \item{object}{
A penalized likelihood model fit obtained by the \code{fit_t_pl} function.
}
  \item{plot}{
Plot of the PC's axes. Default is TRUE (see details).'}

 \item{...}{
Options to be passed through. (e.g., axes=c(1,2), col, pch, cex, mode="cov" or "corr", etc.)}
}


\value{
	a list with the following components
  
  \item{values}{the eigenvalues of the evolutionary variance-covariance matrix}
  \item{scores}{the PC scores}
  \item{loadings}{the component loadings}
  \item{nodes_scores}{the scores for the ancestral states at the nodes (projected on the space of the tips)}
  \item{mean}{the mean/ancestral value used to center the data}
  \item{vectors}{the eigenvectors of the evolutionary variance-covariance matrix}
  }

\details{
\code{phyl.pca_pl} allows computing a phylogenetic principal component analysis (following Revell 2009) using a regularized evolutionary variance-covariance matrix from penalized likelihood models fit to high-dimensional datasets (where the number of variables p is potentially larger than n; see details for the models options in \code{fit_t_pl}). Models estimates are more accurate than maximum likelihood methods, particularly in the high-dimensional case.
Ploting options, the number of axes to display (\code{axes=c(1,2)} is the default), and whether the covariance (\code{mode="cov"}) or correlation (\code{mode="corr"}) should be used can be specified through the ellipsis "..." argument.
}

\note{
	Contrary to conventional PCA, the principal axes of the phylogenetic PCA are not orthogonal, they represent the main axes of (independent) evolutionary changes.}

\references{
Revell, L.J., 2009. Size-correction and principal components for intraspecific comparative studies. Evolution, 63:3258-3268.

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68: 93-116.

}

\author{J. Clavel}

\seealso{
\code{\link{fit_t_pl}},
\code{\link{ancestral}},
\code{\link{GIC.fit_pl.rpanda}}, 
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

# fit a multivariate Pagel lambda model with Penalized likelihood
fit <- fit_t_pl(Y, tree, model="lambda", method="RidgeAlt")

# Perform a phylogenetic PCA using the model fit (Pagel lambda model)
pca_results <- phyl.pca_pl(fit, plot=T) 

# retrieve the scores
head(pca_results$scores)
}
}
