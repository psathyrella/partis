#get Jensen-Shannon divergence	
JSDtree <- function(phylo,meth=c("standard")){
	dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x*log(x/y))
	JSD <- function(x,y) sqrt(0.5*KLD(x,(x+y)/2)+0.5*KLD(y,(x+y)/2))
		matrixColSize <- length(colnames(inMatrix))
		matrixRowSize <- length(rownames(inMatrix))
		colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize) 
	  	inMatrix = apply(inMatrix,1:2,
  		function(x) ifelse (x==0,pseudocount,x))
	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }
#gaussian kernel convolution	 
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    na.omit(x)->x
    if(length(x) == 0)
        stop('too infinite.')
  }
  	kernelG<-function(x, mean=0, sd=1) 
		dnorm(x, mean = mean, sd = sd)
  x<-log(x)		
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

#take square-root of Jensen-Shannon divergence
JSDist <- function(x,y) sqrt(dist.JSD(x,y))
	
#compute eigenvalues for phylogenies and convolve with Gaussian kernel	
	if(meth=="standard"){
		treeNodes <- lapply(phylo,dist.nodes)	 
		treeMats <- lapply(treeNodes,data.matrix)
		treeGraphs <- lapply(treeMats,graph.adjacency,weighted=T)
		treeLaplacian <- lapply(treeGraphs,graph.laplacian,
					normalized=F)			
		treeEigen <- lapply(treeLaplacian,eigen,
				symmetric=TRUE,only.values=TRUE)	
	m<-c()
	d<-c()
		for(i in 1:length(treeEigen)){
			m[[i]]<-subset(treeEigen[[i]]$values,
				treeEigen[[i]]$values>=1)
			d[[i]]<-dens(m[[i]])	
		}
	Ds<-c()
		for(i in 1:length(d)){
			Ds<-as.data.frame(cbind(Ds,d[[i]]$x))
			}
			if (is.null(names(phylo))) {colnames(Ds) <- seq(1,length(d),1)}
			else {colnames(Ds) <- names(phylo)}

			
	JSD<-as.matrix(JSDist(Ds))	
	}	
	
		if(meth=="normal1"){
	treeNodes <- lapply(phylo,dist.nodes)	 
	treeMats <- lapply(treeNodes,data.matrix)
	treeGraphs <- lapply(treeMats,graph.adjacency,weighted=T)
	treeLaplacian <- lapply(treeGraphs,graph.laplacian,
					normalized=T)			
	treeEigen <- lapply(treeLaplacian,eigen,
				symmetric=TRUE,only.values=TRUE)	
	m<-c()
	d<-c()
		for(i in 1:length(treeEigen)){
			m[[i]]<-subset(treeEigen[[i]]$values,
				treeEigen[[i]]$values>=0)
			d[[i]]<-dens(m[[i]])	
		}
	Ds<-c()
		for(i in 1:length(d)){
			Ds<-as.data.frame(cbind(Ds,d[[i]]$x))
			}
			
			if (is.null(names(phylo))) {colnames(Ds) <- seq(1,length(d),1)}
			else {colnames(Ds) <- names(phylo)}
	
	JSD<-as.matrix(JSDist(abs(Ds)))
	}
	
	if(meth=="normal2"){
		treeNodes <- lapply(phylo,dist.nodes)	 
		treeMats <- lapply(treeNodes,data.matrix)
		treeGraphs <- lapply(treeMats,graph.adjacency,weighted=T)
		treeLaplacian <- lapply(treeGraphs,graph.laplacian,
					normalized=F)			
		treeEigen <- lapply(treeLaplacian,eigen,
				symmetric=TRUE,only.values=TRUE)
	m<-c()
	d<-c()	
		for(i in 1:length(treeEigen)){
			m[[i]]<-subset(treeEigen[[i]]$values,
				treeEigen[[i]]$values>=0)
			d[[i]]<-dens(m[[i]]/length(m[[i]]))	
		}
	Ds<-c()
		for(i in 1:length(d)){
			Ds<-as.data.frame(cbind(Ds,d[[i]]$x))
			}
			
			if (is.null(names(phylo))) {colnames(Ds) <- seq(1,length(d),1)}
			else {colnames(Ds) <- names(phylo)}
	
	JSD<-as.matrix(JSDist(abs(Ds)))	
}

#print matrix		

return(JSD)

}	
