\name{likelihood_coal_var}
\alias{likelihood_coal_var}
\title{
Likelihood of a birth-death model using a coalescent approch
}
\description{
Computes the likelihood of a phylogeny under the expanding diversity model with potentially time-varying rates and potentially missing extant species to a phylogeny. Notations follow Morlon et al. PloSB 2010.
}
\usage{
likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0, beta, N0, pos = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Vtimes}{
a vector of branching times (sorted from present to past)
}
  \item{ntips}{
number of species in the phylogeny
}
  \item{lamb0}{
the speciation rate at present
}
  \item{alpha}{
the parameter controlling the exponential variation in speciation rate.  
%Model writes:
%\deqn{\lambda(t) = \lambda_0 e^{\alpha t}}{\lambda(t) = \lambda_0 exp(\alpha t)}
}
  \item{mu0}{
the extinction rate at present
}
  \item{beta}{
the parameter controlling the exponential variation in extinction rate.  
%Model writes:
%\deqn{\mu(t) = \mu_0 e^{\beta t}}{\mu(t) = \mu_0 exp(\beta t)}
}
  \item{N0}{
the number of extanct species
}
  \item{pos}{
logical: should be set to FALSE only to not enforce positive speciation and extinction ratess
}
}

\value{
  a list containing the following components:
  \item{res}{the loglikelihood value of the phylogeny, given the parameters}
  \item{all}{vector of all the individual loglikelihood values corresponding to each branching event}
}

\details{Time runs from the present to the past. Hence, a positive alpha (for example) means that the speciation rate declines from past to present.
%Given a phylogeny, the function computes the likelihood of the distances \eqn{t_1,...,t_{k-2}}{t[1],...,t[k-2]} between
%nodes in that phylogeny.
%If we denote \eqn{t_i}{t[i]} the time-length between nodes i and i+1, and \eqn{v_i}{v[i]} the time-length between node i+1 and the present,
%the likelihood writes as following:
%\deqn{L(t_i) = i(i+1)\frac{\lambda(v_i)}{N(v_i)} exp(-i(i+1) \int_{v_i-t_i}^{v_i} \frac{\lambda(t)}{N(t)})}{L(t[i]) = i(i+1)(lambda(v[i]) / N(v[i])) exp(-i(i+1) integrate(\lambda(t)/N(t) dt, v[i]-t[i], v[i]) )}
%The expression is valid if \eqn{\lambda(t) >= \mu(t)}.
%With fixed speciation and extinction rate functions, log-likelihood is evaluated for each time-length \eqn{t_i}{t[i]} between nodes i and i+1 of the clade, and summed to give the total log-likelihood.

%The likelihood function accounts \eqn{\lambda, \mu}, which might be of different types and may be resumed hereafter:
%\itemize{
%\item model 3: \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0}{\lambda(t) = \lambda_0, \mu(t) = \mu_0}
%(with cst.lamb=TRUE, cst.mu=TRUE)
%\item model 4a: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0}
%(with cst.lamb=FALSE, cst.mu=TRUE)
%\item model 4b: \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0, \mu(t) = \mu_0 exp(\beta t)}
%(with cst.lamb=TRUE, cst.mu=FALSE)
%\item model 4c: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \epsilon \lambda(t)}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \epsilon \lambda(t)}
%(with cst.lamb=FALSE, cst.mu=FALSE, fix.eps=TRUE)
%\item model 4d: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0 exp(\beta t)}
%(with cst.lamb=FALSE, cst.mu=FALSE, fix.eps=FALSE)
%\item model 5: \eqn{\lambda(t) = \lambda_0, \mu(t) = 0}{\lambda(t) = \lambda_0, \mu(t) = 0}
%(with cst.lamb=TRUE, mu.0=TRUE)
%\item model 6: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = 0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = 0}
%(with cst.lamb=FALSE, mu.0=TRUE)
%}
%where \eqn{\lambda_0, \alpha, \mu_0, \beta, \epsilon} are respectively:
%\itemize{
%\item speciation rate at present;
%\item exponential variation in speciation rate;
%\item extinction rate at present;
%\item exponential variation in extinction rate;
%\item extinction fraction (extinction rate/speciation rate).
%}
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the dynamics of diversification: a coalescent approach, PLoS B 8(9): e1000493
}

\author{H Morlon}

\examples{
data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
lamb0 <- 0.1
alpha <- 0.001
mu0<-0
beta<-0
ntips <- Ntip(Cetacea)
N0 <- 89
likelihood <- likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0, beta, N0)
}
