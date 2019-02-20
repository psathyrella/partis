.VCV.rescale.DDexp_geog<-function(phylo,sig2,rate,geo.object,check=FALSE){


paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS


nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
for (i in 1:dim(phylo$edge)[1]){
	nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
}
nodeDist<-c(nodeDist,max(heights))
nodeDiff<-diff(nodeDist)

if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
	node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
	node.order<-node.order+length(phylo$tip.label)
	old.edge<-phylo$edge
	phylo$edge[,1]<-node.order
	for(j in 1:length(phylo$edge[,2])){
		if(phylo$edge[j,2]>length(phylo$tip.label)){
			#match number order in old edge
			#lookup value in new edge
			#replace with value
			phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
			}
		}
	nodeDist<-vector()
	for (i in 1:dim(phylo$edge)[1]){
		nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
		}
	nodeDist<-c(nodeDist,max(heights))
	nodeDiff<-diff(nodeDist)
}

mat<-matrix(nrow=0, ncol=3)
counter_three_letters <- 0
for(i in 1:phylo$Nnode){
	other<-phylo$edge[phylo$edge[,1]==i+length(phylo$tip.label), 2]
	for(b in other){
		int<-matrix(ncol=3)
		int[1]<-i+length(phylo$tip.label)
		if(b>length(phylo$tip.label)){
			counter_three_letters <- counter_three_letters + 1
			int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
			int[3]<-b
			} else {
			int[2]<-phylo$tip.label[[b]]
			int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
			}
		mat<-rbind(mat,int)
		}
	}
		
newDist<-geo.object$times
newDiff<-geo.object$spans
geography.object<-geo.object$geography.object
if(any(nodeDiff==0)){stop("VCV.rescale cannot handle trees with two or more nodes occurring at exactly the same time")}
if(length(geography.object)!=length(newDiff)){stop("The number of sympatry/allopatry matrices does not equal the number of time periods")}
nat<-list()
for(i in 1:length(newDiff)){
	nat[[i]]<-rownames(geography.object[[i]])
	}

if(length(geography.object)>0){
count=1
while(count<=length(nat)){
	if(any(is.na(match(rownames(geography.object[[i]]),unlist(nat[[i]]))))){
		stop("ERROR: Names of geography.object do not match nat list")}
		count=count+1
		}}

V=vcv.phylo(phylo)
N=c(2:(length(newDiff)+1))


if(check==TRUE){
ov<-vector()
for(i in 1:length(phylo$tip.label)){
	for(j in 1:length(phylo$tip.label)){
		sij=V[i,j]
		int<-vector()
		for(m in 1:length(N)){
			if(rownames(V)[i]%in%unlist(nat[[m]])){ #this means that the lineage is present
				geog.int<-sum(geography.object[[m]][match(rownames(V)[i],rownames(geography.object[[m]])),])			
				inta<-(sig2*exp(rate*geog.int))*(max(sij-newDist[m],0)-max(sij-newDist[m+1],0)) 
				int<-c(int,inta)
			} else{
				prev.branch<-mat[mat[,3]==mat[mat[,2]==rownames(V)[i],1],2]
				while(prev.branch%in%unlist(nat[[m]])==FALSE){
					prev.branch<-mat[mat[,3]==mat[mat[,2]==prev.branch,1],2]
					}
				geog.int<-sum(geography.object[[m]][match(prev.branch,rownames(geography.object[[m]])),])			
				inta<-(sig2*exp(rate*geog.int))*(max(sij-newDist[m],0)-max(sij-newDist[m+1],0)) 
				int<-c(int,inta)
			}			
		}
		ov<-c(ov,sum(int))
		}
		}
		
V2<-matrix(ov,nrow=length(phylo$tip.label),byrow=TRUE)
rownames(V2)<-rownames(V)
colnames(V2)<-colnames(V)

if(!isSymmetric(V2)){stop("output VCV not symmetric, can't support current geography matrix")} else{
return(V2)
}
} else {

V2=matrix(nrow=length(phylo$tip.label),ncol=length(phylo$tip.label))
for(i in 1:length(phylo$tip.label)){
	for(j in 1:length(phylo$tip.label)){
		if((lower.tri(V,diag=TRUE)[i,j]==TRUE)){
		sij=V[i,j]
		int<-vector()
		for(m in 1:length(N)){
			if(rownames(V)[i]%in%unlist(nat[[m]])){ #this means that the lineage is present
				geog.int<-sum(geography.object[[m]][match(rownames(V)[i],rownames(geography.object[[m]])),])			
				inta<-(sig2*exp(rate*geog.int))*(max(sij-newDist[m],0)-max(sij-newDist[m+1],0)) 
				int<-c(int,inta)
			} else{
				prev.branch<-mat[mat[,3]==mat[mat[,2]==rownames(V)[i],1],2]
				while(prev.branch%in%unlist(nat[[m]])==FALSE){
					prev.branch<-mat[mat[,3]==mat[mat[,2]==prev.branch,1],2]
					}
				geog.int<-sum(geography.object[[m]][match(prev.branch,rownames(geography.object[[m]])),])			
				inta<-(sig2*exp(rate*geog.int))*(max(sij-newDist[m],0)-max(sij-newDist[m+1],0)) 
				int<-c(int,inta)
			}
		}
		V2[i,j]<-sum(int)				
		}
		}}
		
				
V2[upper.tri(V2)==TRUE]<-t(V2)[upper.tri(t(V2))==TRUE]
rownames(V2)<-rownames(V)
colnames(V2)<-colnames(V)

return(V2)
}		
}