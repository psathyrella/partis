\name{plot_prob_dtt}
\alias{plot_prob_dtt}

\title{
  Plot diversity through time with confidence intervals.
}
\description{
  Plots confidence intervals of the estimated number of species through time using a matrix of probabilities given by the function 'prob_dtt'. 
}
\usage{
  plot_prob_dtt(mat, grain =0.1, plot.prob = TRUE, 
                plot.mean = TRUE, int = TRUE, plot.bound=FALSE,
                conf = 0.95, add = FALSE, col.mean = "red", col.bound = "blue",
                lty="solid", lwd=1)
}

\arguments{
  \item{mat}{
    matrix of probabilities, with species numbers as rows and times as columns with rownames and colnames set to the values of each. 
}
  \item{grain}{
    the upper limit of a range of probabilities plotted in a gray scale (lower limit is zero). Higher probabilities are plotted in black. Default value is 0.1.
}
  \item{plot.prob}{
    logical: set to TRUE (default value) to plot the probabilities.
}
  \item{plot.mean}{
    logical: set to TRUE (default value) to plot a line for the mean.
}
  \item{plot.bound}{
    logical: set to TRUE to plot the bounds of the confidence interval, int must be set to TRUE.
}
  \item{int}{
    logical: set to TRUE (default value) to plot a confidence interval.
}
  \item{conf}{
    confidence level. The default value is 0.95.
}
  \item{add}{
    logical: set to TRUE to add the plot on an existing graph.
}
  \item{col.mean}{
    color of the line for the mean.
}
  \item{col.bound}{
    color of the confidence interval bounds
}
  \item{lty}{
    style of the line for the mean (if added on a current plot)
}
  \item{lwd}{
    the line width, a positive number (default to 1)
}
}
\value{
  Plot representing how the estimated number of species vary through time with confidence intervals. The darker is the plot, the higher is the probability.
}

\details{  
  The function assumes that the matrix of probabilities 'mat' has species numbers as rows and times as columns with rownames and colnames set to the values of each.
  
  'Grain' must be between 0 and 1. If the plot is too pale 'grain' should be diminished (and inversely if the plot is too dark)
}

\author{O.Billaud, T.L.Parsons, D.S.Moen, H.Morlon}

\references{
  Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record. Proc. Nat. Acad. Sci. 108: 16327-16332.
  
  Billaud, O., Moen, D. S., Parsons, T. L., Morlon, H. (under review) Estimating Diversity Through Time using Molecular Phylogenies: Old and Species-Poor Frog Families are the Remnants of a Diverse Past.  
	
}

\seealso{
  \code{\link{fit_bd}}, \code{\link{plot_dtt}}, \code{\link{prob_dtt}}
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

# Plot Diversity through time
plot_prob_dtt(prob)
}