\name{CreateGeobyClassObject}
\alias{CreateGeobyClassObject}

\title{
Create merged biogeography-by-class object
}
\description{

Create a merged biogeography-by-class object to be passed to fit_t_comp_subgroup using a stochastic map created from any model in BioGeoBEARS (see documentation in BioGeoBEARS package) and a simmap object from phytools (see documentation in phytools package).

}
\usage{

CreateGeobyClassObject(phylo,simmap,trim.class,ana.events,clado.events,
	stratified=FALSE,rnd=5)

}

\arguments{
  \item{phylo}{
the object of type 'phylo' (see ape documentation) used to build ancestral range stochastic maps in BioGeoBEARS
}
  \item{simmap}{
a phylo object created using \code{\link{make.simmap}} in phytools
}
  \item{trim.class}{
category in the simmap object that represents the subgroup of interest (see Details and Examples)
}
  \item{ana.events}{
the "ana.events" table produced in BioGeoBEARS that lists anagenetic events in the stochastic map
}
  \item{clado.events}{
the "clado.events" table produced in BioGeoBEARS that lists cladogenetic events in the stochastic map
}
  \item{stratified}{
logical indicating whether the ancestral biogeography stochastic map was built from a stratified analysis in BioGeoBEARS
}
  \item{rnd}{
an integer value indicating the number of decimals to which values should be rounded in order to reconcile class and geo.objects (default is 5)
}



}
\value{

Returns a list with the following components:
\item{map}{a \code{simmap} object with phylogeny trimmed to subgroup of interest (including all branches determined to belong to that subgroup)}

\item{geography.object}{a list with the following components:}
		\item{geography.matrix}{a list of matrices specifying both sympatry & group membership (==1) or allopatry and/or non-membership in the focal subgroup (==0) for each species pair for each internode interval (see Details)}
		\item{times}{a vector containing the time since the root of the tree at which nodes or changes in biogeographyXsubgroup membership occur (used internally in other functions)}
		\item{spans}{a vector specifying the distances between times (used internally in other functions)}
}

\details{

This function merges a class object (which reconstructs group membership through time) and a stochastic map of ancestral biogeography (to reconstruct sympatry through time),
such that lineages can only interact when they belong to the same subgroup AND are sympatric.

This allows fitting models of competition where only sympatric members of a subgroup can compete (e.g., all lineages that share similar diets or habitats).

This function should be used to format the geography object so that it can be correctly
passed to the numerical integration performed in \code{\link{fit_t_comp_subgroup}}. 

}

\references{

	Drury, J., Tobias, J., Burns, K., Mason, N., Shultz, A., and Morlon, H. in review. Contrasting impacts of competition on ecological and social trait evolution in songbirds. \emph{PLOS Biology}.

	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}}

\seealso{
\code{\link{fit_t_comp_subgroup}}, \code{\link{CreateGeoObject_BioGeoBEARS}} , \code{\link{CreateClassObject}} 
}

\examples{

data(BGB.examples)

\donttest{

Canidae.phylo<-BGB.examples$Canidae.phylo
dummy.group<-c(rep("B",3),rep("A",12),rep("B",2),rep("A",6),rep("B",5),rep("A",6))
names(dummy.group)<-Canidae.phylo$tip.label

Canidae.simmap<-phytools::make.simmap(Canidae.phylo,dummy.group)

#build GeobyClass object with "A" as the focal group

Canidae.geobyclass.object<-CreateGeobyClassObject(phylo=Canidae.phylo,simmap=Canidae.simmap, 
trim.class="A",ana.events=BGB.examples$Canidae.ana.events, 
clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE, rnd=5)
	
plotSimmap(Canidae.geobyclass.object$map)
}

}
