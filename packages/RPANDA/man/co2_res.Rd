\name{co2_res}
\alias{co2_res}
\docType{data}
\title{co2 data since the beginning of the Cenozoic}
\description{Atmospheric co2 data since the beginning of the Cenozoic}
\usage{data(co2_res)}
\details{
  Implied co2 data since the beginning of the Cenozoic taken from Hansen et al., (2013). The data are the amount of co2 in ppm reuquired to yield observed global temperature throughout the Cenozoic:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{co2}}{a numeric vector corresponding to the estimated co2 at that age}
  }
}

\source{
Steeman ME et al.(2009) Radiation of extant cetaceans driven by restructuring of the oceans, \emph{Syst Biol} 58:573-585
}

\references{
Steeman ME et al.(2009) Radiation of extant cetaceans driven by restructuring of the oceans \emph{Syst Biol} 58:573-585

Morlon, H., Parsons, T.L., Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332

Condamine, F.L., Rolland, J., Morlon, H. (2013) Macroevolutionary perspectives to environmental change \emph{Eco Lett} 16: 72-85
}
\examples{
data(Cetacea)
print(Cetacea)
plot(Cetacea)
}
\keyword{datasets}
