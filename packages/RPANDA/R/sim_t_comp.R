sim_t_comp<-function(phylo,pars,root.value,Nsegments=1000,model="MC,DDexp,DDlin"){
#return error if non-ultrametric tree
if(phylo$Nnode!=(length(phylo$tip.label)-1)){stop("phylo object must be ultrametric")}
if(length(pars)!=2){stop("pars must be a vector with a value for sig2 and either S for MC model, r for DDexp model, or b for DDlin model")}
if(model=="MC"){sterm=pars[2]}
if(model=="DDlin"){slope=pars[2]}
if(model=="DDexp"){r=pars[2]}
sig2=pars[1]

paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS

ROOT<-root.value

##define a few useful variables
##nodeDist, the distance between each node and the root
##nodeDiff, the distance between each node and the previous node (times for integration)

nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
for (i in 1:dim(phylo$edge)[1]){
	nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
}
nodeDist<-c(nodeDist,max(heights))
nodeDiff<-diff(nodeDist)
###label the branches for each segment of tree to be integrated and identify the node at which the branch terminates

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


##now come up with a list of branches that exist for each time segment (+1 at each branching in this version, which means tree can't have any polytomies or speciation events at exactly the same time)
nat<-list()
for(i in 1:length(nodeDiff)){
	if(i==1){
	nat[[i]]<-list(mat[mat[,1]==(length(phylo$tip.label)+i),2])} else {
	IN<-vector()
	P<-mat[as.numeric(mat[,1])<=(length(phylo$tip.label)+i),c(2,3)]
	IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(length(phylo$tip.label)+i),1])
	nat[[i]]<-list(IN)
	}
	}	
for(i in 2:length(nodeDiff)){			##THIS LOOP checks for an error
	if(length(unlist(nat[[i]]))!=(length(unlist(nat[[i-1]]))+1)){
		print(paste("ERROR at node",i+length(phylo$tip.label)))
		}	
	}

if(is.na(match(model,c("MC","DDexp","DDlin")))){stop("model not specified correctly, must be 'MC','DDexp', or'DDlin'")}
	
if(model=="DDexp"){
	if((sig2/(sig2*exp(r*length(phylo$tip.label))))>=1e3){warning("r parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}
simvalueDDexp<-function(sig2,start.value,seglen,r,branches){
	x<-start.value+rnorm(1,0,sqrt(sig2*exp(r*branches)*seglen))
	return(x)
	}
	}	
		
if(model=="DDlin"){
	if(sig2+(slope*length(phylo$tip.label))<=0){warning("b parameter leads to sig2=0 by the tips, consider changing")}
simvalueDDlin<-function(sig2,start.value,seglen,slope,branches){
	x<-start.value+rnorm(1,0,sqrt(max((sig2+(slope*branches)),0)*seglen))
	return(x)
	}
	}
	
if(model=="MC"){
simvalueMC<-function(sig2,sterm,mu,start.value,seglen){
	x<-start.value+sterm*(mu-start.value)*seglen+rnorm(1,0,sqrt(sig2*seglen))
	return(x)
	}
	}

N<-Nsegments #number of smaller segments to divide a tree into
seglength<-tail(nodeDist,n=1)/N


masterbranch<-list()
masterseg<-list()
for(i in 1:phylo$Nnode){ ##for each node interval
	if(i==1){mu<-ROOT} #initialize mu value at ROOT (change if another root value is desired)
	branchespresent<-length(unlist(nat[[i]]))
	masterbranch[[i]]<-as.list(rep(ROOT,branchespresent))
	names(masterbranch[[i]])<-unlist(nat[[i]])
	if(i>1){
	#update starting values to be the ending values from last time
	for(l in 1:length(masterbranch[[i]])){
		if(!is.na(match(names(masterbranch[[i]][l]),names(masterbranch[[i-1]])))){ #if the name of the element matches branch present in previous iteration
		masterbranch[[i]][l]<-sapply(masterbranch[[i-1]][match(names(masterbranch[[i]][l]),names(masterbranch[[i-1]]))],tail,n=1)		
		} else {
		##look up branch from which descended
		prev.branch<-mat[mat[,3]==mat[mat[,2]==names(masterbranch[[i]][l]),1],2]
		##look up starting value
		masterbranch[[i]][l]<-sapply(masterbranch[[i-1]][match(prev.branch,names(masterbranch[[i-1]]))],tail,n=1)
		}
		}
	}
	masterseg[[i]]<-as.list(rep(0,branchespresent))
	names(masterseg[[i]])<-unlist(nat[[i]])
	for(m in 1:length(masterseg[[i]])){
		if(i==1){
		if(seglength<nodeDiff[i]){
		masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],seq(seglength,nodeDiff[i],seglength))} else{
		masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],nodeDiff[i])	
		}} else{
		masterseg[[i]][[m]]<-seq(nodeDist[i],nodeDist[i+1],seglength)
			}
		if(tail(masterseg[[i]][[m]],n=1)!=nodeDiff[i]){
			masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],nodeDist[i+1])
		}
	}
		
	for(k in 1:length(diff(masterseg[[i]][[1]]))){ ##for each of the segments between intervals
		for(j in 1:branchespresent){## for each of the branches present
			##extract starting value
			start.value<-tail(masterbranch[[i]][[j]],n=1) #replace 1 with matching function
			#run simvalue function with mu (where to get value of mu and update?)
			segsize<-diff(masterseg[[i]][[1]])[k]
			mu.value<-tail(mu,n=1)
			if(model=="DDexp" ){
			sv<-simvalueDDexp(sig2,start.value,segsize,r,branchespresent)
			}
			if(model=="DDlin" ){
			sv<-simvalueDDlin(sig2,start.value,segsize,slope,branchespresent)
			}
			if(model=="MC"){
			sv<-simvalueMC(sig2,sterm,mu.value,start.value,segsize)}
			masterbranch[[i]][[j]]<-c(masterbranch[[i]][[j]],sv)				
			}
		#update mu, add to VECTOR (ie mu=c(mu,UPDATE))
		mu=c(mu,mean(sapply(masterbranch[[i]],tail,n=1)))# the last element in each element of branchlist/branches present ##this assumes that the first element is updated for each i
		}
		}
return(sapply(masterbranch[[phylo$Nnode]],tail,n=1))
}

