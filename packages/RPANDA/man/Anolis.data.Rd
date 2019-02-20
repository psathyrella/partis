\name{Anolis.data}
\alias{Anolis.data}

\title{
Anolis dataset
}
\description{
Phylogeny, trait data, and geography.object for a subclade of Greater Antillean \emph{Anolis} lizards. 
}
\usage{
data(Anolis.data)
}

\details{
	Illustrative phylogeny trimmed from the maximum clade credibility tree of Mahler et al. 2013, corresponding phylogenetic principal component data from Mahler et al. 2013, and biogeography data from Mahler & Ingram 2014 (in the form of a geography object, as detailed in the CreateGeoObject help file). 
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Mahler, D.L., Ingram, T., Revell, L., and Losos, J. 2013. Exceptional convergence on the macroevolutionary landscape in island lizard radiations. \emph{Science}. 341:292-295.
	
	Mahler, D.L. and Ingram, T. 2014. Phylogenetic comparative methods for studying clade-wide convergence. In \emph{Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology}, ed. L. Garamszegi. pp.425-450.
	}
	
\seealso{\code{\link{CreateGeoObject}}}


\examples{
data(Anolis.data)
plot(Anolis.data$phylo)
print(Anolis.data$data)
print(Anolis.data$geography.object)
}
