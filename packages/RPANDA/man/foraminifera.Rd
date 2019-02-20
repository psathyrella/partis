\name{foraminifera}
\alias{foraminifera}
\docType{data}
\title{Foraminifera diversity since the Jurassic}
\description{Foraminifera fossil diversity since the Jurassic}
\usage{data(foraminifera)}
\details{
  Foraminifera fossil diversity since the Jurassic compiled from the Neptune Database (Lazarus, 1994) and Paleobiology Database (https://paleobiodb.org/). Diversity curves are estimated at the genus level using shareholder quorum subsampling (Alroy, 2010) at two-million-year bins. The format is a dataframe with the two following variables:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{foraminifera}}{a numeric vector corresponding to the estimated foraminifera change at that age}
  }
}
\references{
Lazarus, D. (1994) Neptune: A marine micropaleontology database \emph{Mathematical Geology 26:817–832}

Alroy, J. (2010) Geographical, environmental and intrinsic biotic controls on Phanerozoic marine diversification: Controls on phanerozoic marine diversification \emph{Palaeontology 53:1211–1235}

}
\examples{
data(foraminifera)
plot(foraminifera)
}
\keyword{datasets}