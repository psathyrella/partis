##this version allows for a tree to not have sequentially numbered nodes (e.g., trees loaded with read.tree function)
##this version updates the covariance equation to correctly deal with covariance between allopatric lineages
##this version also permits non-island types of biogeography

require(geiger)
require(phytools)
require(deSolve)

.VCV.rescale.geog<-function(phylo,sigma,alpha,sterm,geo.object){
##these are the three parameters that the equations use to produce the variance-covariance matrix, and will ultimately be estimated with ML
parameters<-c(a=alpha,b=sigma,s=sterm) 
if(any(grepl("___",phylo$tip.label))|any(grepl("-",phylo$tip.label))|any(grepl("*",phylo$tip.label))|any(grepl("/",phylo$tip.label))|any(grepl("^",phylo$tip.label))|any(grepl("+",phylo$tip.label))){stop("script will not work with '___', '-', '+', '*','/', or '^' in any tip labels; remove these characters")}
if(!is.binary.tree(phylo)){stop("tree must not contain any polytomies")}
if(sum(phylo$edge.length<0)>0){stop("tree cannot have negative branch lengths")}
if(!is.ultrametric(phylo)){stop("tree must be ultrametric; current verson cannot handle fossil taxa (in development)")}

paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS

root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
totlen<-max(heights)
nodeDist<-c(as.numeric(sort(max(branching.times(phylo))-branching.times(phylo))),totlen)	
nodeDiff<-diff(nodeDist)
if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
old.labels<-as.numeric(names(sort(branching.times(phylo),decreasing=TRUE)))
if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
	checkmat<-cbind(old.labels,seq(root,length(phylo$tip.label)+phylo$Nnode))
	old.edge<-phylo$edge
	for(j in 1:phylo$Nnode){phylo$edge[which(old.edge==checkmat[j,1])]<-checkmat[j,2]}
	}

newDist<-geo.object$times
newDiff<-geo.object$spans
geography.object<-geo.object$geography.object
if(any(nodeDiff==0)){stop("VCV.rescale cannot handle trees with two or more nodes occurring at exactly the same time")}
if(length(geography.object)!=length(newDiff)){stop("The number of sympatry/allopatry matrices does not equal the number of time periods")}

	
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
nat<-list()
for(i in 1:length(newDiff)){
	nat[[i]]<-rownames(geography.object[[i]])
	}
	
#### NOW DEFINE ODEs for each interval and numerically integrate from root to tip, one interval at at a time ###
output<-list() ##initialize storage of results from ODEs 
env.results<-new.env(size = 1, parent = emptyenv())

for(i in 1:length(newDiff)){ ##for each interval between branches...

##var.list is the list of terms for which there is a variance term at each step in time
##cov.list is the list of covariance terms
var.list<-unlist(nat[[i]])
cov.list<-apply(t(combn(var.list,2)),1,paste,collapse="___")
Cmat<-matrix(nrow=length(var.list),ncol=length(var.list)) ##"Cmat" is a lower triangle matrix with each of the variance and covariance terms, used to define equations below
diag(Cmat)<-var.list
Cmat[lower.tri(Cmat)]<-cov.list

##define ODEs for deSOLVE
	len <- length(var.list)
	dim <- len * (len + 1) / 2
	ind <- matrix(0, nrow = len, ncol=len, byrow=TRUE)
	indI <- vector(mode = "integer", length = dim)
	indJ <- vector(mode = "integer", length = dim)
	for(k in 1:len) {
	    ind[[k,k]] <- k
	    indI[[k]] <- k
	    indJ[[k]] <- k
	}
	counter <- len
	for(l in 1:(len-1)) {
	    for(k in (l+1):len) {
	        counter <- counter + 1
	        ind[[k,l]] <- counter
	        ind[[l,k]] <- counter
	        indI[[counter]] <- k
	        indJ[[counter]] <- l
	    }
	}
	coefAlpha <- (2*(len-1)/len * parameters['s'] + 2 * parameters['a'])
	coefBeta <- parameters['s'] / len
	rhs <- vector(mode = "numeric", length = dim)
	for(m in 1:len) {
	    rhs[m] <- parameters['b']
	}
	A <- diag(x = -coefAlpha, nrow = dim, ncol=dim)
	if(parameters['s']!=0){
	updated.alpha<-vector()	
	for(term in 1:length(indJ)){
		jj<-indJ[term]
		ii<-indI[term]
		nj=sum(geography.object[[i]][,jj])
		ni=sum(geography.object[[i]][,ii])
		int<-((1-(1/(2*nj))-(1/(2*ni)))/((len-1)/len)*(coefAlpha-2* parameters['a']))+2* parameters['a']
		updated.alpha<-c(updated.alpha, int)
		}	
	diag(A)<--updated.alpha
	}
	coefBetaV <- rep(coefBeta, times = len)
	for(term in 1:len){
		int<-coefBeta*(len/sum(geography.object[[i]][,term]))
		coefBetaV[term]<-int
		}	
	for(ind_jk in 1:dim) {
		    j <- indI[[ind_jk]]
		    k <- indJ[[ind_jk]]
		    geog<-geography.object[[i]][match(ind[k,],ind)]
			coefB<-rep(coefBetaV[k],len) #is the right coefBetaV term being repped
	    	coefB[k]<-0
		    A[ind_jk, ind[j, ] ] <- (A[ind_jk, ind[j, ] ] + coefB)*geog 
	    	geog<-geography.object[[i]][match(ind[j,],ind)]
    		coefB<-rep(coefBetaV[j],len)
    		coefB[j]<-0
	    	A[ind_jk, ind[k, ] ] <- (A[ind_jk, ind[k, ] ] + coefB)*geog
			}
	diag(A)<-sapply(diag(A),function(x) {if(x==0){x<--(2* parameters['a'])}else{x<-x}})
	
 # return the rate of change
 sttrt=paste("d",c(var.list,cov.list),sep="",collapse=",")
 sttat=paste("list(c(",sttrt,"))",sep="")

######	STARTING VALUES  ######
#This next block of code provides the starting values for the ODEs by identifying the proper value from the previous integration
if(i==1){ #if it is the first iteration, the starting values for all terms are 0
start=c(paste(c(var.list,cov.list),"=0",sep="",collapse=","))
state.new<-paste("c(",start,")",sep="")
state<-eval(parse(text=state.new))
} else {
#initialize starting value named vector
start=c(paste(c(var.list,cov.list),"=0",sep="",collapse=","))
state.new<-paste("c(",start,")",sep="")
state<-eval(parse(text=state.new))
termlist<-c(var.list,cov.list)
for(l in 1:length(termlist)){
term<-termlist[l]
if(exists(term, envir=env.results)){ #if the term to be numerically integrated was present in the previous generation(branch) of the tree
       state[l]<-get(term, envir=env.results)
	}  
else { #if term is not present in previous generation/branch of the tree, lookup which branch it descends from
	if(l <= length(var.list)){ ##this loop looks up values for variance terms using the "mat" matrix that is produced earlier
		prev.branch<-mat[mat[,3]==mat[mat[,2]==term,1],2]
		state[l]<-get(prev.branch, envir=env.results)
	}
	else{ ##for covariance terms that weren't present in previous generation/branch, it is necessary to decompose covariance term, look up terms appropriately, and re-assemble them in the correct order so they match the lower triangle format used in definition of terms
		firstterm<-strsplit(term,"___")[[1]][1]
		scondterm<-strsplit(term,"___")[[1]][2]
		rearr<-paste(text=scondterm,"___",firstterm,sep="")
		if(exists(rearr,envir=env.results)){
			state[l]<-get(rearr, envir=env.results)
		}else{
		if(exists(firstterm, envir=env.results) || exists(scondterm, envir=env.results)){
			if(exists(firstterm, envir=env.results)){
				prev.branch<-mat[mat[,3]==mat[mat[,2]==scondterm,1],2]
				prev.branch2<-firstterm} 
			else{  ##NOTE prev.branch2 just identifies which term was present previously
				prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
				prev.branch2<-scondterm
				}
		prev.term<-paste(prev.branch2,prev.branch,sep="___")
		if(! exists(prev.term, envir=env.results)){
			prev.term<-paste(prev.branch,prev.branch2,sep="___")
			}
		state[l]<-get(prev.term, envir=env.results)
		} 				
	else {
		##neither first nor second term was present in the previous generation/branch, so the starting value is a variance term from t-1
		prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
		state[l]<-get(prev.branch, envir=env.results)
	}
	}
	}
	}
	}
	}

ou <- function(t, state, parameters) {
  dX <- A %*% state + rhs
  return (list(c(dX)))
}

##NOW, run numerical integration for given time step and append to a list
output<-ode(y=state,times=c(0,newDiff[i]),func=ou,parms=NULL)
colN <- colnames(output)
env.results<-new.env(size = length(colN), parent = emptyenv())
for (k in 2:length(colN)){
  assign(colN[[k]], output[[2,k]], envir = env.results)
}
}

Vou<-matrix(nrow=length(phylo$tip.label),ncol=length(phylo$tip.label))
rownames(Vou)<-unlist(nat[[length(newDiff)]])
colnames(Vou)<-unlist(nat[[length(newDiff)]])
diag(Vou)<-output[2,][2:(length(phylo$tip.label)+1)]
for(j in (length(phylo$tip.label)+2):length(output[2,])){

	string<-strsplit(names(output[2,j]),"___")
	string<-unlist(string)
	one<-match(string[1],colnames(Vou))
	two<-match(string[2],colnames(Vou))
	Vou[one,two]<-output[2,j]
	Vou[two,one]<-output[2,j]
}

return(Vou)
}

