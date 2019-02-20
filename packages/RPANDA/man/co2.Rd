\name{co2}
\alias{co2}
\docType{data}
\title{co2 data since the Jurassic}
\description{Atmospheric co2 data since the Jurassic}
\usage{data(co2)}
\details{
  Atmospheric co2 data since the Jurassic taken from Mayhew et al., (2008, 2012) and derived from the GeoCarb-III model (Berner and Kothavala, 2001). The data are eported as the ratio of the mass of co2 at time t to that at present. The format is a dataframe with the two following variables:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{co2}}{a numeric vector corresponding to the estimated co2 at that age}
  }
}
\references{
Mayhew, P.J., Jenkins, G.B., Benton, T.G. (2008) A long-term association between global temperature and biodiversity, origination and extinction in the fossil record \emph{Proceedings of the Royal Society B 275:47-53}

Mayhew, P.J., Bell, M.A., Benton, T.G, McGowan, A.J. (2012) Biodiversity tracks temperature over time \emph{109:15141-15145}

Berner R.A., Kothavala, Z. (2001) GEOCARB III: A revised model of atmospheric CO2 over
Phanerozoic time \emph{Am J Sci 301:182–204}

}
\examples{
data(co2)
plot(co2)
}
\keyword{datasets}