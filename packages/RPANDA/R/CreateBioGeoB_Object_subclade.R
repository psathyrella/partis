.CreateBioGeoB_Object_subclade<-function(anc.phylo,subclade.phylo,ana.events,clado.events,nat.only=FALSE){
	phylo<-subclade.phylo
	subclade.tips<-phylo$tip.label
	paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
	paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	totlen<-length(phylo$tip.label)
	root <-totlen  + 1
	heights<-nodeHeights(phylo)
	totheight=max(nodeHeights(phylo))
	totheight.anc=max(nodeHeights(anc.phylo))
	subclade.root=getMRCA(anc.phylo,subclade.tips)
	tips = clado.events$node[which(clado.events$label%in%subclade.tips)]
	nn<-unique(apply(combn(tips,2),2,function(x)getMRCA(anc.phylo,x)))
	M<-mrca(anc.phylo,full=T)
	desc=unique(as.vector(M[which(rownames(M)%in%c(tips,nn)),]))
	desc=desc[c(which(desc>=subclade.root),which(desc<=length(anc.phylo$tip.label)))]
	subnodes=c(desc[which(desc>length(anc.phylo$tip.label))])
	ghost.nodes=desc[which(!desc%in%c(nn,tips))]
	ana.events<-ana.events[which(ana.events$node%in%desc[which(desc!=subclade.root)]),]
	clado.events<-clado.events[which(clado.events$node%in%subnodes),]
	nodes=unique(round(clado.events[,9],8))
	events=unique(round(totheight.anc-ana.events$abs_event_time,8))
	nodeDist=sort(c(nodes,events,totheight.anc))
	nodeDiff=diff(nodeDist)
	old.phylo<-phylo
	old.edge<-anc.phylo$edge[which(anc.phylo$edge[,1]%in%nn),]
	other.edge<-phylo$edge
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	old.labels<-as.numeric(names(sort(branching.times(phylo),decreasing=TRUE)))
	if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
		checkmat<-cbind(old.labels,seq(root,length(phylo$tip.label)+phylo$Nnode))
		for(j in 1:phylo$Nnode){phylo$edge[which(other.edge==checkmat[j,1])]<-checkmat[j,2]}
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
	desc.mat<-matrix(nrow=length(desc), ncol=2)
	desc.mat[,1]<-desc
	for(i in 1:length(desc)){
		if(desc.mat[i,1]%in%ghost.nodes){
			flag=0
			j=1
			while(flag==0){
			term=anc.phylo$edge[which(anc.phylo$edge[,1]==desc.mat[i,1]),2][j]
			if(term%in%desc){
			flag=1
			}else{
			j=2
			}}
			if(term%in%tips){
				desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
			} 
			if(term%in%nn){
				desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
			}
			term2=term 
			while(term2%in%ghost.nodes){
				flag=0
				j=1
				while(flag==0){
				term=anc.phylo$edge[which(anc.phylo$edge[,1]==term2),2][j]
				if(term%in%desc){
				flag=1
				}else{
				j=2
				}}
				if(term%in%tips){
					desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
				} 
				if(term%in%nn){
					desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
				} 
				term2=term
				}
			} 
		if(desc.mat[i,1]%in%tips){
			desc.mat[i,2]<-which(anc.phylo$tip.label[desc.mat[i,1]]==phylo$tip.label)
		}
		if(desc.mat[i,1]%in%nn){
				desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==desc.mat[i,1]),1][1]
		}
	}
	
	
	#while(any(desc.mat[,2]%in%ghost.nodes)){
	#	nxt<-desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]
	#	desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]<-desc.mat[match(nxt ,desc.mat[,1]),2]
	#}
	
	nat<-list()
	nodecount=1
	for(i in 1:length(nodeDiff)){
		if(nodeDist[i]%in%nodes){ #a new species appears, this is a node
			if(clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]%in%nn){## node has two descendants because of subtree trimming
				IN<-vector()
				node<-totlen+nodecount
				P<-mat[as.numeric(mat[,1])<=(node),c(2,3)]
				IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(node),1])
				#need to find a way to look up old node
				#oldnode=old.edge[which(phylo$edge[,1]==node)][1]
				delbr<-clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]
				left=clado.events[which(clado.events[,1]==delbr),15]
				left=desc.mat[which(desc.mat[,1]==left),2]
				right=clado.events[which(clado.events[,1]==delbr),16]
				right=desc.mat[which(desc.mat[,1]==right),2]
				m<-regexpr("->",clado.events[which(clado.events[,1]==delbr),20])
				regs<-regmatches(clado.events[which(clado.events[,1]==delbr),20],m,invert=TRUE)[[1]][2]
				p<-regexpr(",",regs)	
				sides<-regmatches(regs,p,invert=TRUE)
				lside<-sides[[1]][1]
				rside<-sides[[1]][2]
				if(i==1){
					geo.vector<-vector(length=length(IN))
					if(left<=totlen){
						geo.vector[which(IN==phylo$tip.label[left])]<-lside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
					}
					if(right<=totlen){
						geo.vector[which(IN==phylo$tip.label[right])]<-rside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
					}
					nat[[i]]<-cbind(IN,geo.vector)
				} else {
					ogv<-geo.vector			##need to update old geo.vector
					geo.vector<-vector(length=length(IN))
					terms<-which(nat[[i-1]][,1]%in%IN) #this should give the item numbers of old things
					geo.vector[match(nat[[i-1]][,1][terms],IN)]<-ogv[terms]
					if(left<=totlen){
						geo.vector[which(IN==phylo$tip.label[left])]<-lside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
					}
					if(right<=totlen){
						geo.vector[which(IN==phylo$tip.label[right])]<-rside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
					}
					nat[[i]]<-cbind(IN,geo.vector)			
				}} else{
				IN<-nat[[i-1]][,1]
				delbr<-clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]
				left=clado.events[which(clado.events[,1]==delbr),15]
				right=clado.events[which(clado.events[,1]==delbr),16]
				m<-regexpr("->",clado.events[which(clado.events[,1]==delbr),20])
				regs<-regmatches(clado.events[which(clado.events[,1]==delbr),20],m,invert=TRUE)[[1]][2]
				p<-regexpr(",",regs)	
				sides<-regmatches(regs,p,invert=TRUE)
				lside<-sides[[1]][1]
				rside<-sides[[1]][2]
				if(left%in%desc){ #left branch should be recorded
					hold<-desc.mat[which(desc.mat[,1]==left),2]
					if(hold<=totlen){
						elno<-which(IN==phylo$tip.label[hold])
					}else{
						elno<-which(IN==mat[,2][mat[,3]==hold])
					}
					geo.vector[elno]<-lside
				} else{#right branch should be recorded
					hold<-desc.mat[which(desc.mat[,1]==right),2]
					if(hold<=totlen){
						elno<-which(IN==phylo$tip.label[hold])
					}else{
						elno<-which(IN==mat[,2][mat[,3]==hold])
					}
					geo.vector[elno]<-rside
				}
				nat[[i]]<-cbind(IN,geo.vector)			
		}} else{ # this is a branch change
				IN<-nat[[i-1]][,1]
				#identify which branch changes
				##note should eventually move this up so subtraction doesn't happen every time(for speed up)
				delbr<-which(nodeDist[i]==round(totheight.anc-ana.events$abs_event_time,8))
				if(length(delbr)>1){stop("more than one events at same time")}
				hold<-desc.mat[which(desc.mat[,1]==ana.events[delbr,]$nodenum_at_top_of_branch),2]
				if(hold<=totlen){
					elno<-which(IN==phylo$tip.label[hold])
				}else{
					elno<-which(IN==mat[,2][mat[,3]==hold])
				}
				#identify new state
				geo.vector[elno]<-ana.events[delbr,]$new_rangetxt
				nat[[i]]<-cbind(IN,geo.vector)			
			}
	if(nodeDist[i+1]%in%nodes && clado.events[which(round(clado.events[,9],8)==round(nodeDist[i+1],8)),1]%in%nn){nodecount=nodecount+1}
	if(any(nat[[i]][,2]==FALSE)){stop("ERROR:'FALSE' recorded as range")}
	}			
	
	if(nat.only==TRUE){
		return(nat[[length(nat)]])
	} else{
	
	coll.vector<-vector()
	count=1
	geography.matrix<-list()
		for(i in 1:length(nodeDiff)){ 
			timei<-unlist(nat[[i]])
			var.list<-timei[,1]
			geo.list<-strsplit(timei[,2],"")
			len=length(var.list)
			int.mat<-matrix(nrow=len,ncol=len)
			rownames(int.mat)<-var.list
			colnames(int.mat)<-var.list
			diag(int.mat)<-1
			for(j in 1:length(var.list)){
				for(k in 1:length(var.list)){
					if(j>k){
						int.mat[j,k]<-ifelse(any(geo.list[[j]]%in%geo.list[[k]]),1,0)
					}
				}
			}
			int.mat[upper.tri(int.mat)==TRUE]<-t(int.mat)[upper.tri(t(int.mat))==TRUE]
			if(count>1 && !is.character(all.equal(int.mat,geography.matrix[[count-1]]))){coll.vector<-c(coll.vector,i)} else{	
			geography.matrix[[count]]<-int.mat
			count=count+1
			}
		}
	nodeDist<-nodeDist[which(!nodeDist%in%nodeDist[coll.vector])]
	nodeDist<-nodeDist-min(nodeDist)
	nodeDiff<-diff(nodeDist)	
return(list(geography.object=geography.matrix,times=nodeDist,spans=nodeDiff))
}
}
