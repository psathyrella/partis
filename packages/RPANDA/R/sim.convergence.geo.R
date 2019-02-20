
sim.convergence.geo<-function(phylo,pars, Nsegments=2500, plot=FALSE, geo.object){
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
  root <- length(phylo$tip.label) + 1
  heights<-phytools::nodeHeights(phylo)
  totlen<-max(heights)
  len<-root-1
  nodeDist<-c(as.numeric(sort(max(ape::branching.times(phylo))-ape::branching.times(phylo))),totlen)	
  nodeDiff<-diff(nodeDist)
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  old.labels<-as.numeric(names(sort(ape::branching.times(phylo),decreasing=TRUE)))
  old.edge<-phylo$edge
  if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
  	checkmat<-cbind(old.labels,seq(root,len+phylo$Nnode))
  	old.edge<-phylo$edge
  	for(j in 1:phylo$Nnode){phylo$edge[which(old.edge==checkmat[j,1])]<-checkmat[j,2]}
  	}
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+len, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+len
      if(b>len){
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
  branchesPresent <- rep(NA, length(nodeDiff)) 
  for(i in 1:length(nodeDiff)){
    if(i==1){
        nat[[i]]<-mat[mat[,1]==(len+i),2]
      } else {
        IN<-vector()
        P<-mat[as.numeric(mat[,1])<=(len+i),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(len+i),1])
        nat[[i]]<-IN
      }
    branchesPresent[i] = length(nat[[i]])
  }	
    
  masterbranch<-list()
  masterbranch2<-list()

  segmentSize <- rep(NA, phylo$Nnode)
  mappings <- list()
  segsize = sum(nodeDiff)/Nsegments
  newDist<-geo.object$times
  newDiff<-geo.object$spans
  geography.object<-geo.object$geography.object

  for(i in 1:phylo$Nnode){ ##for each node interval
    
    nati<-nat[[i]]
    mappings[[i]] = rep(NA, branchesPresent[i])
    
    if(i==1){
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
    }else{      
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
      for (j in 1:branchesPresent[i]){
        if(length(which(nati[j] == pastNati)) > 0){ 
          mappings[[i]][j] = which(nati[j] == pastNati)
        }else{
          mappings[[i]][j] = which(mat[mat[,3]==mat[mat[,2]==nati[j],1],2] == pastNati)
        }
      }     
    }   
    pastNati <- nati
  }
  
  ## Looping over the parameters
  
  out <- matrix( nrow = nrow(pars), ncol = len)
  out2 <- matrix( nrow = nrow(pars), ncol = len)
  
  for (p in 1:nrow(pars)){
    sig2 = pars[p,1]
    max = pars[p,2]
    alpha = pars[p,3]
    root.value = pars[p,4]
    psi=pars[p,5]
    theta=pars[p,6]
    
    
    timecount=1      
    for(i in 1:phylo$Nnode){   
        traitMat <- matrix(nrow = branchesPresent[i], ncol = segmentSize[i]+1)
        traitMat2 <- matrix(nrow = branchesPresent[i], ncol = segmentSize[i]+1)
        
        if (i == 1){
          traitMat[,1] = root.value
          traitMat2[,1] = root.value
        }else{
          traitMat[,1] = masterbranch[[i-1]][mappings[[i]],(segmentSize[i-1]+1)] #added +1 here since segmentSize[i-1] is penultimate column, not the last one
		  traitMat2[,1] = masterbranch2[[i-1]][mappings[[i]],(segmentSize[i-1]+1)] 
        
        }
        
        tempInd<- 1:branchesPresent[i] # hack to have fast selection of not k, seemed to be faster than a call to which()
        
        for(k in 1:segmentSize[i]){        
          	for(j in 1:branchesPresent[i]){
	          	elms<-which(geography.object[[timecount]][match(unlist(nat[[i]])[j],rownames(geography.object[[timecount]])),]==1)##these are elements that are equal to 1 (sympatric lineages);if length(elms)=0 the repulsion component just falls out
	            diffMeOthers <- traitMat2[j,k] - traitMat2[elms[which(elms!=j)],k]
	            sign1<-sign(traitMat[j,k] - traitMat[elms[which(elms!=j)],k]) ##Alternative: sign1<-ifelse(diffMeOthers>0,1,-1), but diff seems to be faster
	            temp1 = max *segsize  #old version:temp1 = 1/(length(elms)) * max *segsize #
	        	if(abs(temp1)==Inf){temp1=0}
      		    temp2 = sqrt(sig2*segsize)    
	        	traitMat[j,k+1]<-traitMat[j,k] + temp1 * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2)
	            traitMat2[j,k+1]<-traitMat2[j,k] + psi*(theta-traitMat2[j,k])*segsize + rnorm(1,0,temp2)
              	
				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),8)>=round(newDist[timecount+1],8),timecount+1,timecount)}
					###loop for last segment size (to preserve exact branch lengths)
					if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
						segsizeT= nodeDiff[i]%%segsize
	            		diffMeOthers <- traitMat2[j,k] - traitMat2[elms[which(elms!=j)],k]
						temp1B =  max *segsizeT # 1/length(elms) * max *segsizeT  # #old version 
	        			if(abs(temp1B)==Inf){temp1B=0}
	        			temp2B = sqrt(sig2*segsizeT)
	            		sign1<-sign(traitMat[j,k] - traitMat[elms[which(elms!=j)],k]) # Alternative: sign1<-ifelse(diffMeOthers>0,1,-1), but diff seems to be faster
	            		traitMat[j,k+1]<-traitMat[j,k]  + temp1B * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2B)
	            		traitMat2[j,k+1]<-traitMat2[j,k] + psi*(theta-traitMat2[j,k])*segsizeT + rnorm(1,0,temp2B)
					 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),8)>=round(newDist[timecount+1],8),timecount+1,timecount)}
						}           
	        	}
    	}
    	
        masterbranch[[i]] = traitMat
        masterbranch2[[i]] = traitMat2
      }
    out[p,] = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])      
    if(p==1){colnames(out)<-unlist(nat[[i]])}
    
    out2[p,] = as.vector(masterbranch2[[phylo$Nnode]][, tail(segmentSize,n=1)+1])      
    if(p==1){colnames(out2)<-unlist(nat[[i]])}			
			
  }

  if(plot==TRUE){
    print("plotting last simulated dataset")
    par(mfrow=c(2,1))
    M=seq(0,sum(nodeDiff),length=sum(segmentSize))
    O=list()
    O2=list()
    for(i in 1:length(nodeDiff)){
    O[[i]]<-data.frame(seq(nodeDist[[i]],nodeDist[[i+1]],length=(segmentSize[[i]]+1)),as.data.frame(t(masterbranch[[i]])))
    O2[[i]]<-data.frame(seq(nodeDist[[i]],nodeDist[[i+1]],length=(segmentSize[[i]]+1)),as.data.frame(t(masterbranch2[[i]])))
    }
    t.plot<-plot(M,1:length(M),col="white", ylim=c(range(sapply(masterbranch,range))), xlab="time", ylab="Value")
    for(i in 1:length(nodeDiff)){
    	for(j in 1:(ncol(O[[i]])-1)){
    		lines(O[[i]][,1],O[[i]][,j+1])
    	}	
    }
    t2.plot<-plot(M,1:length(M),col="white", ylim=c(range(sapply(masterbranch2,range))), xlab="time", ylab="Value")
    for(i in 1:length(nodeDiff)){
    	for(j in 1:(ncol(O2[[i]])-1)){
    		lines(O2[[i]][,1],O2[[i]][,j+1])
    	}	
    }
    
  }

  return(list(trait1=out,non.focal=out2))
}



