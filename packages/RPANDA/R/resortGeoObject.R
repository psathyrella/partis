.resortGeoObject<-function(phylo,geo.object){
	gmat<-geo.object$geography.object
	if(any(grepl("___",phylo$tip.label))|any(grepl("-",phylo$tip.label))|any(grepl("/",phylo$tip.label))){stop("script will not work with '___', '-', '+', '*','/', or '^' in any tip labels; remove these characters")}
	paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
	paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
	totlen<-length(phylo$tip.label)
	root <-totlen  + 1
	heights<-nodeHeights(phylo)
	for (i in 1:dim(phylo$edge)[1]){
		nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
	}
	nodeDist<-c(nodeDist,max(heights))
	nodeDiff<-diff(nodeDist)
	flag=0
	if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
		node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
		node.order<-node.order+totlen
		old.edge<-phylo$edge
		old.phylo<-phylo
		phylo$edge[,1]<-node.order
		for(j in 1:length(phylo$edge[,2])){
			if(phylo$edge[j,2]>totlen){
				#match number order in old edge
				#lookup value in new edge
				#replace with value
				phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
				}
			}
		nodeDist<-vector()
		for (i in 1:dim(phylo$edge)[1]){
			nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
			}
		nodeDist<-c(nodeDist,max(heights))
		nodeDiff<-diff(nodeDist)
		flag=1
	}
	
	
	tips<-1:length(phylo$tip.label) 
	order.vec<-rep(0,length=dim(phylo$edge)[1])

	counter=1
	for(i in (totlen+1):(totlen+phylo$Nnode)){
	 	for(j in 1:2){
	 		if(order.vec[which(phylo$edge[,1]==i)][j]==0){
	  			order.vec[which(phylo$edge[,1]==i)][j]<-counter
				m<-phylo$edge[which(phylo$edge[,1]==i)[j],2]
				while(!m%in%tips){
	 				order.vec[which(phylo$edge[,1]==m)][1]<-counter
	 				m<-phylo$edge[which(phylo$edge[,1]==m)[1],2]
	 			 }
	 			counter=counter+1 
	 			} 
	 			}
	 		}
		
	
	newedge<-cbind(phylo$edge,order.vec)
	
	mat<-matrix(nrow=0, ncol=4)
	counter_three_letters <- 0
	for(i in 1:phylo$Nnode){
		other<-newedge[newedge[,1]==i+totlen, 2]
		for(b in other){
			int<-matrix(ncol=4)
			int[1]<-i+totlen
			if(b>totlen){
				counter_three_letters <- counter_three_letters + 1
				int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
				int[3]<-b
				int[4]<-newedge[which(newedge[,2]==b),3]
				} else {
				int[2]<-phylo$tip.label[b]
				int[3]<-b
				int[4]<-newedge[which(newedge[,2]==b),3] 
				}
			mat<-rbind(mat,int)
			}
		}
	
	sorted.gmat<-list()
	for(i in 1:length(gmat)){
		sorted.gmat[[i]]<-gmat[[i]][order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4])),order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4]))]
	}
	
	
	return(list(geography.object=sorted.gmat,times=geo.object$times,spans=geo.object$spans))
}