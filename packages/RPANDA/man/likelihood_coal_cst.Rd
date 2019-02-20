\name{likelihood_coal_cst}
\alias{likelihood_coal_cst}

\title{
Likelihood of a phylogeny under the equilibrium diversity model
}
\description{
Computes the likelihood of a phylogeny under the equilibrium diversity model with potentially time-varying rates and potentially missing extant species. Notations follow Morlon et al. PloSB 2010.
}
\usage{
likelihood_coal_cst(Vtimes, ntips, tau0, gamma, N0)
}
\arguments{
  \item{Vtimes}{
a vector of branching times (sorted from present to past)
}
  \item{ntips}{
the number of tips in the phylogeny
}
  \item{tau0}{
the turnover rate at present
}
  \item{gamma}{
the parameter controlling the exponential variation in turnover rate. With gamma=0, the turnover rate is constant over time.
}
  \item{N0}{
the number of extant species
}
}

\value{
  a list containing the following components:
  \item{res}{the loglikelihood value of the phylogeny, given tau0 and gamma}
  \item{all}{vector of all the individual loglikelihood values corresponding to each branching event}
}

\details{Time runs from the present to the past. Hence, a positive gamma (for example) means that the turnover rate declines from past to present.
%Given a phylogeny, the function computes the likelihood of the distances \eqn{t_1,...,t_{k-2}} between
%nodes in that phylogeny using a constant diversity hypothesis.
%If we denote \eqn{t_i} the time-length between nodes i and i+1, \eqn{v_i} is the time-length between node i+1 and the present,
%the likelihood writes as following:
%
%\eqn{L(t_i) = i(i+1) \frac{\tau_0}{N_0} exp(-i(i+1) \frac{\tau_0}{N_0 \gamma} exp(\gamma v_i) (1-exp(-\gamma t_i)) ), \gamma \neq 0}
%
%\eqn{L(t_i) = i(i+1) \frac{\tau_0}{N_0} exp(-i(i+1) \frac{\tau_0}{N_0} t_i), \gamma = 0}
%
%where:
%\itemize{
%\item \eqn{\tau_0} is the turnover rate at present time;
%\item \eqn{\gamma} is exponential variation in turnover rate.
%}
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the dynamics of diversification: a coalescent approach, PLoS B 8(9): e1000493
}

\author{H Morlon}

\examples{
data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
tau0 <- 0.1
gamma <- 0.001
ntips <- Ntip(Cetacea)
N0 <- 89
likelihood <- likelihood_coal_cst(Vtimes,ntips,tau0,gamma,N0)
}
