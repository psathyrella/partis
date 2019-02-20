\name{plot_fit_bd}
\alias{plot_fit_bd}

\title{
Plot speciation, extinction & net diversification rate functions of a fitted model
}
\description{
Plot estimated speciation, extinction & net diversification rates through time
}
\usage{
plot_fit_bd(fit.bd, tot_time)
}

\arguments{
  \item{fit.bd}{
an object of class 'fit.bd', output of the 'fit_bd' function
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
}


\value{
Plots representing how the estimated speciation, extinction & net diversification rate functions vary through time
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
result <- fit_bd(Balaenopteridae,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     expo.lamb = TRUE, fix.mu=TRUE)
# plot fitted rates
#plot_fit_bd(result, tot_time)
}
