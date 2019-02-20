CreateGeoObject<-function(phylo,map){
  	if(any(grepl("___",phylo$tip.label))){stop("script will not work with '___' in tip labels; remove extra underscores")}
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
	mat<-matrix(nrow=0, ncol=3)
	counter_three_letters <- 0
	for(i in 1:phylo$Nnode){
		other<-phylo$edge[phylo$edge[,1]==i+totlen, 2]
		for(b in other){
			int<-matrix(ncol=3)
			int[1]<-i+totlen
			if(b>totlen){
				counter_three_letters <- counter_three_letters + 1
				int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
				int[3]<-b
				} else {
				int[2]<-phylo$tip.label[[b]]
				int[3]<-0 
				}
			mat<-rbind(mat,int)
			}
		}
####				
####
	if(any(class(map)=="phylo")){
	for(i in 1:length(mat[,1])){
		if(mat[i,3]==0){
		mat[i,3]<-as.character(match(mat[i,2],phylo$tip.label))
	}}	
	maps.object<-map$maps
	brchange<-vapply(maps.object,function(x) any(length(names(x))>1),1)
	if(sum(brchange)!=0){
	newtimes<-vector()
	for(i in 1:length(brchange)){
		if(brchange[i]!=0){
			#look up branch in stochastic map
			intlen<-length(maps.object[[i]])
			inttime<-vector()
			for(j in 1:(intlen-1)){
				inttime<-c(inttime,maps.object[[i]][j])
				nt<-as.numeric(heights[i,1]+sum(inttime))
				newtimes<-c(newtimes,nt)
				}
		}
	}
	old.Dist<-nodeDist
	old.Diff<-nodeDiff
	nodeDist<-sort(c(nodeDist,newtimes))
	nodeDiff<-diff(nodeDist)
	}
	nat<-list()
	nodecount=1
	for(i in 1:length(nodeDiff)){
		if(i==1){
			hold.m<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
			int<-dim(hold.m)[1]
			for(m in 1:int){
				hold.m[m,2]<-names(maps.object[[which(phylo$edge[,2]==as.numeric(hold.m[m,2]))]])[1]
			}
			nat[[i]]<-hold.m
			} else {
			if(nodeDist[i]%in%old.Dist){ #a new species appears, this is a node
				P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
				hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
				int<-dim(hold.m)[1]
				for(m in 1:int){
					iden<-which(phylo$edge[,2]==as.numeric(hold.m[m,2]))
					if(brchange[iden]==0 || !(hold.m[m,1]%in%nat[[i-1]][,1])){
						hold.m[m,2]<-names(maps.object[[iden]])[1]
					} else{
						num=1
						while(round(nodeDist[i+1]-old.Dist[phylo$edge[iden,1]-totlen],6)>round(sum(maps.object[[iden]][1:num]),6)){
						num=num+1}
						hold.m[m,2]<-names(maps.object[[iden]])[num]
					}
				}
				nat[[i]]<-hold.m
			} else{ # this is a branch change
				P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
				hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
				int<-dim(hold.m)[1]
				for(m in 1:int){
					iden<-which(phylo$edge[,2]==as.numeric(hold.m[m,2]))
					if(brchange[iden]==0 | !(hold.m[m,1]%in%nat[[i-1]][,1])){
						hold.m[m,2]<-names(maps.object[[iden]])[1]
					} else{					
						num=1
						while(round(nodeDist[i+1]-old.Dist[phylo$edge[iden,1]-totlen],6)>round(sum(maps.object[[iden]][1:num]),6)){
						num=num+1}
						hold.m[m,2]<-names(maps.object[[iden]])[num]
					}
				}
			nat[[i]]<-hold.m
			}
		}
	if(nodeDist[i+1]%in%old.Dist){nodecount=nodecount+1}
	}	
	
	geography.matrix<-list()
		for(i in 1:length(nodeDiff)){ 
			timei<-unlist(nat[[i]])
			var.list<-timei[,1]
			len=length(var.list)
			int.mat<-matrix(nrow=len,ncol=len)
			rownames(int.mat)<-var.list
			colnames(int.mat)<-var.list
			diag(int.mat)<-1
			for(j in 1:length(var.list)){
				for(k in 1:length(var.list)){
					if((lower.tri(int.mat,diag=TRUE)[j,k]==TRUE)){
						int.mat[j,k]<-ifelse(timei[match(var.list[j],timei[,1]),2]==timei[match(var.list[k],timei[,1]),2],1,0)
					}
				}
			}
			int.mat[upper.tri(int.mat)==TRUE]<-t(int.mat)[upper.tri(t(int.mat))==TRUE]
			geography.matrix[[i]]<-int.mat
		}}
####		
####		
		if(is.matrix(map)){
			nat<-list()
			for(i in 1:length(nodeDiff)){
				if(i==1){
				nat[[i]]<-list(mat[mat[,1]==(totlen+i),2])} else {
				IN<-vector()
				P<-mat[as.numeric(mat[,1])<=(totlen+i),c(2,3)]
				IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(totlen+i),1])
				nat[[i]]<-list(IN)
				}
				}	
				for(i in 1:length(mat[,1])){
						if(mat[i,3]==0){
						mat[i,3]<-as.character(match(mat[i,2],phylo$tip.label))
					}}
			for(i in 1:length(mat[,1])){
				if(mat[i,3]==0){
				mat[i,3]<-as.character(match(mat[i,2],phylo$tip.label))
			}}	
			
			if(flag==0){
				map.v<-map[match(mat[,3],map[,2]),3]
				}
			if(flag==1){  
				map2<-cbind(phylo$edge,map[,3])
				map.v<-map2[match(mat[,3],map2[,2]),3]
			}
			mat<-data.frame(mat,map.v)
			#colnames(mat)<-c("node.number","descendant.branch","next.node","region")
			geography.matrix<-list()
			for(i in 1:phylo$Nnode){ 
				var.list<-unlist(nat[[i]])
				len=length(var.list)
				int.mat<-matrix(nrow=len,ncol=len)
				rownames(int.mat)<-var.list
				colnames(int.mat)<-var.list
				diag(int.mat)<-1
				for(j in 1:length(var.list)){
					for(k in 1:length(var.list)){
						if((lower.tri(int.mat,diag=TRUE)[j,k]==TRUE)){
							int.mat[j,k]<-ifelse(mat[match(var.list[j],mat[,2]),4]==mat[match(var.list[k],mat[,2]),4],1,0)
						}
					}
				}
				int.mat[upper.tri(int.mat)==TRUE]<-t(int.mat)[upper.tri(t(int.mat))==TRUE]
				geography.matrix[[i]]<-int.mat
			}
		}						
		
		
####
####								
		return(list(geography.object=geography.matrix,times=nodeDist,spans=nodeDiff))
	}
	
