\name{plot_fit_env}
\alias{plot_fit_env}

\title{
Plot speciation, extinction & net diversification rate functions of a fitted environmental model
}
\description{
Plot estimated speciation, extinction & net diversification rates as a function of the environmental data and time
}
\usage{
plot_fit_env(fit.env, env_data, tot_time)
}

\arguments{
  \item{fit.env}{
an object of class 'fit.env', output of the 'fit_env' function
}
  \item{env_data}{
environmental data, given as a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
}


\value{
Plots representing how the estimated speciation, extinction & net diversification rate functions vary as a function of the environmental data & time
}

\author{H Morlon and FL Condamine}

\seealso{
\code{\link{fit_env}}
}

\examples{

data(Balaenopteridae)
tot_time<-max(node.age(Balaenopteridae)$ages)
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate
# with time
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
#result <- fit_env(Balaenopteridae,InfTemp,tot_time,f.lamb,f.mu,
#      lamb_par,mu_par,f=1, fix.mu=TRUE, df=dof, dt=1e-3)

# plot fitted rates
#plot_fit_env(result, InfTemp, tot_time)
}
