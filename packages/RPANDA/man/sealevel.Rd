\name{sealevel}
\alias{sealevel}
\docType{data}
\title{Sea level data since the Jurassic}
\description{Global sea level change since the Jurassic}
\usage{data(sealevel)}
\details{
  Eustatic sea level change since the Jurassic calculated by Miller et al. (2005) from satellite measurements, tide gauges, shoreline markers, reefs, atolls, oxygen isotopes,, the flooding history of continental margins, cratons. The format is a dataframe with the two following variables:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{sea level}}{a numeric vector corresponding to the estimated sea level change at that age}
  }
}
\references{
Miller, K.G., Kominz, M.A., Browning, J.V., Wright, J.D., Mountain, G.S., Katz, M.E., Sugarman, P.J., Cramer, B.S., Christie-Blick, N., Pekar, S.F. (2005) The Phanerozoic Record of Global Sea-Level Change \emph{Science 310:1293-1298}

}
\examples{
data(sealevel)
plot(sealevel)
}
\keyword{datasets}