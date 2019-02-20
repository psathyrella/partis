\name{createModelCoevolution}
\alias{createModelCoevolution}
\title{
Creation of a PhenotypicGMM
}
\description{
Creates an object of class PhenotypicGMM, a subclass of the class PhenotypicModel intended to represent the Generalist Matching Mutualism model of trait evolution on two specific trees.
}
\usage{
createModelCoevolution(tree1, tree2, keyword)
}

\arguments{
\item{tree1}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{tree2}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{keyword}{
a string object. Defaut value "GMM" returns an object of class PhenotypicGMM, which takes advantage of faster distribution computation. Otherwise, a "PhenotypicModel" is returned, and the computation of the tip distribution will take much longer.
}
}

\value{
an object of class "PhenotypicModel" or "PhenotypicGMM".
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

\examples{
#Loading example trees
newick1 <- "(((A:1,B:1):3,(C:3,D:3):1):2,E:6);"
tree1 <- read.tree(text=newick1)
newick2 <- "((X:1.5,Y:1.5):3,Z:4.5);"
tree2 <- read.tree(text=newick2)

#Creating the model
modelGMM <- createModelCoevolution(tree1, tree2)

#Printing basic or full informations on the model definitions
show(modelGMM)
print(modelGMM)

#Simulates tip trait data
dataGMM <- simulateTipData(modelGMM, c(0,0,5,-5, 1, 1), method=2)
}
