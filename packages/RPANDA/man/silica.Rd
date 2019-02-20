\name{silica}
\alias{silica}
\docType{data}
\title{Silica data across the Cenozoic}
\description{Silica weathering ratio across the Cenozoic}
\usage{data(silica)}
\details{
  Silica weathering ratio across the Cenozoic calculated by Cermeno et al. (2015) using  the lithium isotope record of seawater from Misra and Froelich (2012). The format is a dataframe with the two following variables:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{silica weathering ratio}}{a numeric vector corresponding to the estimated CO2 at that age}
  }
}
\references{
 Misra, S., Froelich, P.N. (2012) Lithium isotope history of Cenozoic seawater: Changes in
silicate weathering and reverse weathering. \emph{Science 335(6070):818–823}

Cermeno, P., Falkowski, P.G., Romero, O.E., Schaller, M.F., Vallina, S.M. (2015) Continental erosion and the Cenozoic rise of marine diatoms \emph{Proceedings of the National Academy of Sciences 112:4239-244}

}
\examples{
data(silica)
plot(silica)
}
\keyword{datasets}