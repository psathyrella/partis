\name{InfTemp}
\alias{InfTemp}
\docType{data}
\title{Paleotemperature data across the Cenozoic}
\description{Paleotemperature data across the Cenozoic inferred from delta O18 measurements}
\usage{data(InfTemp)}
\details{
  Paleotemperature data inferred from delta 018 measurements using the equation of Epstein et al. (1953). The format is a dataframe with the two following variables:
  \describe{
    \item{\code{Age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{Temperature}}{a numeric vector corresponding to the inferred temperature at that age}
  }
}
\references{
Epstein, S., Buchsbaum, R., Lowenstam, H.A., Urey, H.C. (1953) Revised carbonate-water isotopic temperature scale \emph{Geol. Soc. Am. Bull. 64: 1315-1326}

Zachos, J.C., Dickens, G.R., Zeebe, R.E. (2008) An early Cenozoic perspective on greenhouse warming and carbon-cycle dynamics \emph{Nature} 451: 279-283 

Condamine, F.L., Rolland, J., Morlon, H. (2013) Macroevolutionary perspectives to environmental change \emph{Eco Lett} 16: 72-85
}
\examples{
data(InfTemp)
plot(InfTemp)
}
\keyword{datasets}
