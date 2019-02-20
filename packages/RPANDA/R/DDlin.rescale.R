.vcv.rescale.DDlin<-function(phylo,sig2,rate){

bt=as.numeric(branching.times(phylo))
nodeDist<-c(bt[1]-sort(bt,decreasing=TRUE),bt[1])
nodeDiff<-diff(nodeDist)
V=vcv.phylo(phylo)
N=c(2:(length(nodeDiff)+1))

vect.unique <- function(ll)  {
  ll[!duplicated(ll)]
}

sij.list<-unique(round(vect.unique(V),digits=6))
sij<-matrix(nrow=length(sij.list),ncol=2)
sij[,1]<-sij.list
for(m in 1:length(sij.list)){
	int<-vector()
	for(i in 1:length(N)){
		inta<-(sig2+(rate*N[i]))*(max(sij[m,1]-nodeDist[i],0)-max(sij[m,1]-nodeDist[i+1],0)) 
		int<-c(int,inta)
		}
	sij[m,2]<-sum(int)
	}
V2<-matrix(sij[,2][match(sapply(V,round,digits=6),round(sij[,1],digits=6))],nrow=nrow(V))
rownames(V2)<-rownames(V)
colnames(V2)<-colnames(V)

V2
}