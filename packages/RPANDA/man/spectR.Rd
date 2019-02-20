\name{spectR}
\alias{spectR}

\title{
Spectral density plot of a phylogeny
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree, identifies the spectral gap, then convolves the eigenvalues with a Gaussian kernel, and plots them alongside all eigenvalues ranked in descending order.
}

\usage{
spectR(phylo, method=c("standard"))
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{method}{
the method used to compute the spectral density, which can either be "standard" or "normal". If set to "standard", computes the unnormalized version of the spectral density. If set to "normal", computes the spectral density normalized to the degree matrix (see the associated paper for an explanation)}
}

\value{
	a list with the following components:
	\item{eigenvalues}{
the vector of eigenvalues} 
	\item{principal_eigenvalue}{
the largest (or principal) eigenvalue of the spectral density profile}
\item{asymmetry}{
the skewness of the spectral density profile} 
	\item{peak_height}{
the largest y-axis value of the spectral density profile}
	\item{eigengap}{
the position of the largest difference between eigenvalues, giving the number of modalities in the tree}
	}
	
\details{
	Note that the eigengap should in principle be computed with the "standard" option}


\references{
	Lewitus, E., Morlon, H., Characterizing and comparing phylogenies from their Laplacian spectrum, bioRxiv doi: http://dx.doi.org/10.1101/026476}

\author{E Lewitus}

\seealso{
\code{\link{plot_spectR}}, \code{\link{JSDtree}}, \code{\link{BICompare}}
}


\examples{
data(Cetacea)
spectR(Cetacea,method="standard")
}
