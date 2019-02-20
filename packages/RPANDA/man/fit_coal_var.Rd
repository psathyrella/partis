\name{fit_coal_var}
\alias{fit_coal_var}
\title{
Fit birth-death model using a coalescent approch
}
\description{
Fits the expanding diversity model with potentially time-varying rates and potentially missing extant species to a phylogeny, by maximum likelihood. 
The implementation allows only exponential time variation of the speciation and extinction rates, although this could be modified using expressions in Morlon et al. PloSB 2010. Notations follow Morlon et al. PLoSB 2010.
}

\usage{
fit_coal_var(phylo, lamb0 = 0.1, alpha = 1, mu0 = 0.01, beta = 0,
             meth = "Nelder-Mead", N0 = 0, cst.lamb = FALSE, cst.mu = FALSE,
             fix.eps = FALSE, mu.0 = FALSE, pos = TRUE)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{lamb0}{
initial value of the speciation rate at present (used by the optimization algorithm)
}
  \item{alpha}{
 initial value of the parameter controlling the exponential variation in speciation rate (used by the optimization algorithm)
%the exponential variation in speciation rate. Model writes \eqn{\lambda(t) = lamb0 e^{\alpha*t}}{lambda(t) = lamb0 exp(alpha t)}
}
  \item{mu0}{
initial value of the extinction rate at present (used by the optimization algorithm)
}
  \item{beta}{
initial value of the parameter controlling the exponential variation in extinction rate.
%Model writes \eqn{\mu(t) = mu0 e^{\beta*t}}{mu(t) = mu0 exp(beta t)}
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
  \item{N0}{
Number of extant species.
With default value(0), N0 is set to the number of tips in the phylogeny. That is, the phylogeny is assumed to be 100\% complete. 
}
  \item{cst.lamb}{
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time, models 3, 4b & 5 in Morlon et al. PloSB 2010) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time, models 3 & 4a in Morlon et al. PloSB 2010) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{fix.eps}{
logical: should be set to TRUE only if the extinction fraction is constant (i.e. does not depend on time, model 4c in Morlon et al. PloSB 2010)
}
  \item{mu.0}{
logical: should be set to TRUE to force the extinction rate to 0 (models 5 & 6 in Morlon et al. PloSB 2010)
}
  \item{pos}{
logical: should be set to FALSE only to not enforce positive speciation and extinction rates
}
}

\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{model.parameters}{the estimated parameters}
}

\details{
The function fits models 3 to 6 from the PloSB 2010 paper. Likelihoods arising from these models are computed using the coalescent approximation and are directly comparable to likelihoods from the fit_coal_cst function, thus allowing to test support for equilibrium versus expanding diversity scenarios.

%The default parameters allow to fit the most general model where speciation and extinction rates vary over time without a fixed extinction fraction and write as exponential models.
%Depending on the model, these parameters include a combination of the speciation rate at present (lamb0), the exponential variation in speciation rate (alpha), the extinction rate at present (mu0), the exponential variation in% extinction rate (beta) and the extinction fraction (extinction rate/speciation rate, eps). See notations in the PLoS Biology 2010 paper.

These models can be fitted using the options specified below:
\itemize{
\item model 3: 
%\eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0}{\lambda(t) = \lambda_0, \mu(t) = \mu_0}
with cst.lamb=TRUE & cst.mu=TRUE
\item model 4a: 
%\eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0}
with cst.lamb=FALSE & cst.mu=TRUE
\item model 4b: 
%\eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0, \mu(t) = \mu_0 exp(\beta t)}
with cst.lamb=TRUE & cst.mu=FALSE
\item model 4c: 
%\eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \epsilon \lambda(t)}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \epsilon \lambda(t)}
with cst.lamb=FALSE, cst.mu=FALSE & fix.eps=TRUE
\item model 4d: 
%\eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0 exp(\beta t)}
with cst.lamb=FALSE, cst.mu=FALSE & fix.eps=FALSE
\item model 5: 
%\eqn{\lambda(t) = \lambda_0, \mu(t) = 0}{\lambda(t) = \lambda_0, \mu(t) = 0}
with cst.lamb=TRUE & mu.0=TRUE
\item model 6: 
%\eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = 0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = 0}
with cst.lamb=FALSE & mu.0=TRUE
}
%where \eqn{\lambda_0, \alpha, \mu_0, \beta, \epsilon} are respectively:
%\itemize{
%\item the speciation rate at present;
%\item the exponential variation in speciation rate;
%\item the extinction rate at present;
%\item the exponential variation in extinction rate;
%\item the extinction fraction (extinction rate/speciation rate).
%}
Time runs from the present to the past. Hence, if alpha is estimated to be positive (for example), this means that the speciation rate decreases from past to present. 
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the dynamics of diversification: a coalescent approach, PLoS B 8(9): e1000493

Morlon, H., Kemps, B., Plotkin, J.B., Brisson, D. (2012) Explosive radiation of a bacterial species group, Evolution, 66: 2577-2586

Morlon, H. (2014) Phylogenetic approaches for studying diversification, Eco Lett, 17:508-525
}

\author{H Morlon}

\seealso{
\code{\link{likelihood_coal_var}}, \code{\link{fit_coal_cst}}
}

\examples{
data(Cetacea)
result <- fit_coal_var(Cetacea, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=89)
print(result)
}
