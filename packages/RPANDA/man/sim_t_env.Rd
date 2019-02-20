\name{sim_t_env}
\alias{sim_t_env}

\title{
Recursive simulation (root-to-tip) of the environmental model 
}
\description{
Simulates datasets for a given phylogeny under the environmental model (see ?fit_t_env)
}
\usage{

sim_t_env(phylo, param, env_data, model, root.value=0, step=0.001, plot=FALSE, ...)

}
\arguments{
  \item{phylo}{
An object of class 'phylo' (see ape documentation)
}
  \item{param}{
A numeric vector of parameters for the user-defined climatic model. For the \emph{EnvExp} and \emph{EnvLin}, there is only two parameters. The first is sigma and the second beta.
}
  \item{env_data}{
Environmental data, given as a time continuous function (see, e.g. splinefun) or a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}

  \item{model}{
The model describing the functional form of variation of the evolutionary rate \eqn{\sigma^2} with time and the environmental variable. Default models are "EnvExp" and "EnvLin" (see details). An user defined function of any functional form can be used (forward in time).
This function has three arguments: the first argument is time; the second argument is the environmental variable; the third argument is a numeric vector of the parameters controlling the time and environmental variation (to be estimated). See the example below.
}
  \item{root.value}{
A number specifying the trait value for the ancestor
}
  \item{step}{
This argument describe the length of the segments to simulate across for the phylogeny. The smaller is the segment, the greater is the accuracy of the simulation at the expense of the computation time. 
}
  \item{plot}{
If TRUE, the simulated process is plotted.
}

  \item{...}{
Arguments to be passed through. For instance, "col" for plot=TRUE. 
}
}


\value{
A named vector with simulated trait values for \eqn{n} species in the phylogeny
}

\details{
The users defined function is simulated forward in time i.e.: from the root to the tips. The speed of the simulations might depend on the value used for the "step" argument. 
It's possible to estimate the traits with the MLE from another fitted object (see the example below).
}

\references{
Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Science,  114(16): 4183-4188.

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{fit_t_env}}  
\code{\link{plot.fit_t.env}}, 
\code{\link{likelihood_t_env}}
}

\examples{

\donttest{
data(Cetacea)
data(InfTemp)


set.seed(123)
# define the parameters
param <- c(0.1, -0.5)
# define the environmental function
my_fun <- function(t, env, param){ param[1]*exp(param[2]*env(t))}

# simulate the trait
trait <- sim_t_env(Cetacea, param=param, env_data=InfTemp, model=my_fun, root.value=0,
                    step=0.001, plot=TRUE)

# fit the model to the simulated trait.
fit <- fit_t_env(Cetacea, trait, env_data=InfTemp, model=my_fun, param=c(0.1,0))
fit

# Then use the results from the previous fit to simulate a new dataset
trait2 <- sim_t_env(Cetacea, param=fit, step=0.001, plot=TRUE)
fit2 <- fit_t_env(Cetacea, trait2, env_data=InfTemp, model=my_fun, param=c(0.1,0))
fit2

# When providing the environmental function:
require(pspline)
spline_result <- sm.spline(x=InfTemp[,1],y=InfTemp[,2], df=50)
env_func <- function(t){predict(spline_result,t)}
t<-unique(InfTemp[,1])

# We build the interpolated smoothing spline function
env_data<-splinefun(t,env_func(t))

# provide the environmental function to simulate the traits
trait3 <- sim_t_env(Cetacea, param=param, env_data=env_data, model=my_fun,
                     root.value=0, step=0.001, plot=TRUE)
fit3 <- fit_t_env(Cetacea, trait3, env_data=InfTemp, model=my_fun, param=c(0.1,0))
fit3

  }
}
