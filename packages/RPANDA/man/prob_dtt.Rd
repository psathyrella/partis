\name{prob_dtt}
\alias{prob_dtt}

\title{
  Confidence intervals of diversity through time
}
\description{
  Returns a matrix of probabilities to have 'm' species at a given time 't' with 'n' observed extant species (complete sampling or not) and 's' species at the root of the phylogeny (s=1 if the tree has a stem, otherwise s=2)
}
\usage{
  prob_dtt(fit.bd, tot_time, time, N0, l=N0, f = l/N0, 
            m = seq(N0), method="simple", lin = FALSE,
           prec = 1000, type = "stem",logged = TRUE)
}

\arguments{
  \item{fit.bd}{
    an object of class 'fit.bd', output of the 'fit_bd' function.
}
  \item{tot_time}{
    the age of the underlying phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{time}{
    vector of times on which the function calculates the probabilities of 'm' species. The function goes forward in time, so that \eqn{t = 0} is the time of the most recent common ancestor.
}
  \item{N0}{
    number of extant species. If all extant species are represented in the phylogeny, N0 is given by length(phylo$tip.label).
}
  \item{l}{
    number of extant species sampled. Default value is N0 (complete sampling).
}
  \item{f}{
    the fraction of extant species included in the phylogeny, given by l/N0.
}
  \item{m}{
    a vector of integers for which we want to know the probability of each value. 
}
  \item{method}{
    reflects which way of computing is choosen. A 'simple' one (quicker) is used when the number of extant species (N0) is known exactly or when the whole phylogeny is sampled (f==1). A 'hard one', much longer, is used when N0 is not known without doubt and f<1. The default value is "simple"" (the other possibility is "hard")
  }
  \item{lin}{
    logical: set to TRUE if \eqn{\lambda} & \eqn{\mu} are fitted with a linear model.
}
  \item{prec}{
    precision (number of bits used) of the computation. The default value is 1000.
}
  \item{type}{
    reflects whether the clade has a stem or not. Options are the default "stem"" and the alternative "crown", which means the tree starts with two species at time 0.
  }
  \item{logged}{
    logical: set to TRUE to log probabilities and factorials as much as possible (required, except perhaps for very small, young clades).
  }
}

\value{
  Matrix of probabilities to have 'm' species at a given time 't' with 'n' observed extant species (complete sampling or not).
}

\details{
  If the sampling fraction is not equal to 1, the function computes with very high numbers. To be sufficiently accurate, the package 'Rmpfr' is used and "prec" is the precision of the computation. Hence, the calculation may take a lot of time. In case of wrong probabilities (negatives or higher than 1 for instance) you should increase the precision.
  
  If the sampling fraction is equal to 1, the function doesn't need the package 'Rmpfr' and simply uses the log of probabilities and factorials (argument "logged"). Thus, computation is faster. 
  
  The matrix columns names go backward in time.
}

\author{O.Billaud, T.L.Parsons, D.S.Moen, H.Morlon}

\references{
  Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record. Proc. Nat. Acad. Sci. 108: 16327-16332.
  
  Billaud, O., Moen, D. S., Parsons, T. L., Morlon, H. (under review) Estimating Diversity Through Time using Molecular Phylogenies: Old and Species-Poor Frog Families are the Remnants of a Diverse Past. 
  
}

\seealso{
  \code{\link{fit_bd}}, \code{\link{plot_dtt}}, \code{\link{plot_prob_dtt}}
}

\examples{
data(Balaenopteridae)
tot_time<-max(node.age(Balaenopteridae)$ages)

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate
# with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.08, 0.01)
mu_par<-c()
result <- fit_bd(Balaenopteridae,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=1,
                     expo.lamb = TRUE, fix.mu=TRUE)
                     
# Compute the matrix of probabilities                     
prob <- prob_dtt(result, tot_time, 1:tot_time, N0=9, type="crown")

# Check that the sums of probabilities are equal to 1
colSums(prob)
}