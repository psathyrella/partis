\name{fit_coal_cst}
\alias{fit_coal_cst}

\title{
Maximum likelihood fit of the equilibrium model
}
\description{
Fits the equilibrium diversity model with potentially time-varying turnover rate and potentially missing extant species to a phylogeny, by maximum likelihood. 
The implementation allows only exponential time variation of the turnover rate, although this could be modified using expressions in Morlon et al. PloSB 2010. Notations follow Morlon et al. PLoSB 2010.
}
\usage{
fit_coal_cst(phylo, tau0 = 1e-2, gamma = 1, cst.rate = FALSE,
             meth = "Nelder-Mead", N0 = 0)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{tau0}{
initial value of the turnover rate at present (used by the optimization algorithm)
}
  \item{gamma}{
initial value of the parameter controlling the exponential variation in turnover rate (used by the optimization algorithm)
}
  \item{cst.rate}{
logical: should be set to TRUE to fit an equilibrium diversity model with time-constant turnover rate (know as the Hey model, model 1 in Morlon et al. PloSB 2010). By default, a model with expontential time-varying rate exponential is fitted (model 2 in Morlon et al. PloSB 2010).
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
  \item{N0}{Number of extant species.
With default value(0), N0 is set to the number of tips in the phylogeny. That is, the phylogeny is assumed to be 100\% complete. 
}
}

\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{tau0}{the estimated turnover rate at present}
  \item{gamma}{the estimated parameter controlling the exponential variation in turnover rate (if cst.rate is FALSE)}
}

\details{
This function fits models 1 (when cst.rate=TRUE) and 2 (when cst.rate=FALSE) from the PloSB 2010 paper. Likelihoods arising from these models are directly comparable to likelihoods from the fit_coal_var function, thus allowing to test support for equilibrium versus expanding diversity scenarios. Time runs from the present to the past. Hence, if gamma is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
%The function fits the constant diversity/exponentially turnover rate model to a given phylogeny, by maximizing the associated likelihood function.
%Assuming we know the total number of species at present, we estimate the turnover rate parameters. 
%In addition to the turnover rate at present/the exponential variation in turnover rate, the log-likelihood and the second order Akaike's Information Criterion are returned.
%The underlying model (Hey/Moran Model) may be resumed hereafter:
%\eqn{ \tau(t) = \tau_0 e^{-\gamma t}}{ \tau(t) = \tau_0 exp(-\gamma t)}
%where \eqn{\tau(t)} is the turnover rate function.
%Under these assumptions, the model is said to be satured diversity with positive extinction.
%For more details, see notations in the PloSB 2010 paper (model 2).
%The meth argument fixes the optimization method. Depending on this choice, one should take
%care of the default values of tau0 (1e-2) and gamma (1). Indeed, the likelihood function has several local maxima.
%To make sure to not be stuck in a local optimum, several inital values might be investigated.
}

\references{
Hey, J. (1992) Using phylogenetic trees to study speciation and extinction, Evolution, 46: 627-640

Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the dynamics of diversification: a coalescent approach, PLoS B, 8(9): e1000493

Morlon, H., Kemps, B., Plotkin, J.B., Brisson, D. (2012) Explosive radiation of a bacterial species group, Evolution, 66: 2577-2586

Morlon, H. (2014) Phylogenetic approaches for studying diversification, Eco Lett, 17:508-525
}

\author{H Morlon}

\seealso{
\code{\link{likelihood_coal_cst}}, \code{\link{fit_coal_var}}
}

\examples{
data(Cetacea)
result <- fit_coal_cst(Cetacea, tau0=1.e-3, gamma=-1, cst.rate=FALSE, N0=89)
print(result)
}
