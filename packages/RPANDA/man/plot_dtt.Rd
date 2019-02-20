\name{plot_dtt}
\alias{plot_dtt}

\title{
Plot diversity through time
}
\description{
Plot the estimated number of species through time
}
\usage{
plot_dtt(fit.bd, tot_time, N0)
}

\arguments{
  \item{fit.bd}{
an object of class 'fit.bd', output of the 'fit_bd' function
}
  \item{tot_time}{
the age of the underlying phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{N0}{
number of extant species. If all extant species are represented in the phylogeny, N0 is given by length(phylo$tip.label)
}

}


\value{
Plot representing how the estimated number of species vary through time
}

\references{
	Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	
	Morlon, H. (2014) Phylogenetic approaches for studying diversification, Eco Lett 17:508-525
}

\author{H Morlon}

\seealso{
\code{\link{fit_bd}}
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

# plot estimated number of species through time
# plot_dtt(result, tot_time, N0=9)
}
