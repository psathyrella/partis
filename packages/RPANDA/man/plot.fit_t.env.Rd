\name{plot.fit_t.env}
\alias{plot.fit_t.env}

\title{
Plot the phenotypic evolutionary rate through time estimated by the \emph{fit_t_env} function
}
\description{
Plot estimated evolutionary rate as a function of the environmental data and time.
}
\usage{

\method{plot}{fit_t.env}(x, steps = 100, ...)

}
\arguments{
  \item{x}{
an object of class 'fit_t.env' obtained from a \emph{fit_t_env} fit.
}
  \item{steps}{
the number of steps from the root to the present used to compute the evolutionary rate \eqn{\sigma2} through time.
}
  \item{...}{
further arguments to be passed to \code{plot}. See ?\code{plot}.
}

}


\value{
	\code{plot.fit_t.env} returns invisibly a list with the following components used in the current plot:
  
  \item{time_steps}{the times steps where the climatic function was evaluated to compute the rate. The number of steps is controlled through the argument \code{steps}.}
  \item{rates}{the estimated evolutionary rate through time estimated at each \code{time_steps}}
  
}


\note{
	All the graphical parameters (see \code{\link{par}}) can be passed through (e.g. line type: \code{lty}, line width: \code{lwd}, color: \code{col} ...)
}
\references{
Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Science, 114(16): 4183-4188.

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{lines.fit_t.env}}, \code{\link{likelihood_t_env}}
}

\examples{
data(Cetacea)
data(InfTemp)

# Simulate a trait with temperature dependence on the Cetacean tree
set.seed(123)
\donttest{

trait <- sim_t_env(Cetacea, param=c(0.1,0.2), env_data=InfTemp, model="EnvExp", 
					root.value=0, step=0.01, plot=TRUE)


## Fit the Environmental-exponential model

result1=fit_t_env(Cetacea, trait, env_data=InfTemp, scale=TRUE)
plot(result1)

# further options
plot(result1, lty=2, lwd=2, col="red")
}

}
