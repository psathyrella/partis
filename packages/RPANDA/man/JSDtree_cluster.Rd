\name{JSDtree_cluster}
\alias{JSDtree_cluster}

\title{
Clustering of phylogenies
}
\description{Clusters phylogenies using hierarchical and k-medoids clustering}
\usage{
JSDtree_cluster(JSDtree,alpha=0.9,draw=T)
}

\arguments{
  \item{JSDtree}{
a matrix of distances between phylogenie pairs, typically the output of the JSDtree function when the distance is measured as the Jensen-Shannon distance 
}
 \item{alpha}{
   	the confidence value for demarcating clusters in the hierarchical clustering plot; the default is 0.9
	}
\item{draw}{
	plot heatmap and hierarchical cluster in new windows
}	
	}

\value{
	plots a heatmap and a hierarchical cluster with bootstrap support, and outputs results of the k-medoids clustering in the form of a list with the following components
	\item{clusters}{the optimal number of clusters around medoids (see pamk documentation)}
	\item{cluster_assignments}{assignments of trees to clusters} 
	\item{cluster_support}{a list with the following components: widths: a table specifying the cluster to which each tree belongs, the neighbor (i.e. most similar) cluster, and the silhouette width of the observation (see silhouette documentation); clus.avg.widths: average silhouette width for each cluster; vg.width: average silhouette width across all clusters} 
		}

\note{
	The k-medoids clustering may not work with fewer than 10 trees}

\author{E Lewitus}

\references{
	Lewitus, E., Morlon, H., Characterizing and comparing phylogenies from their Laplacian spectrum, bioRxiv doi: http://dx.doi.org/10.1101/026476}

\seealso{
\code{\link{JSDtree}}
}

\examples{
trees<-TESS::tess.sim.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
res<-JSDtree(trees)
#JSDtree_cluster(res,alpha=0.9,draw=T)
}
