\name{d13c}
\alias{d13c}
\docType{data}
\title{d13c data since the Jurassic}
\description{Benthic d13c weathering ratio since the Jurassic}
\usage{data(d13c)}
\details{
  Ratio of stable carbon isotopes since the Jurassic calculated by Hannisdal and Peters (2011) and Lazarus et al. (2014) from marine carbonates. The format is a dataframe with the two following variables:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{d13c}}{a numeric vector corresponding to the estimated d13c at that age}
  }
}
\references{
Hannisdal, B., Peters, S.E. (2011) hanerozoic Earth system evolution and marine biodiversity \emph{Science 334:1121-1124}

Lazarus, D., Barron, J., Renaudie, J., Diver, P., Turke, A. (2014) Cenozoic Planktonic Marine Diatom Diversity and Correlation to Climate Change \emph{PLoS ONE 9:e84857}

}
\examples{
data(d13c)
plot(d13c)
}
\keyword{datasets}
