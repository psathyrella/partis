#This is a script that trims a stochastic map to preserve only those branches containing a given category
#The inputs are:
#map (simmap object that is a stochastic map of discrete traits used to identify guild of competitors) 
#trim.class (discrete category identifying species in 'data' (e.g., herbivores))

.trimSimmap<-function(map,trim.class){
		trc=trim.class
		smap<-map
		
		repeat{
		
		tot.len<-length(smap$tip.label)
		node = tot.len + 1
		tbdropped<-c()
		
		while(node <= (tot.len+smap$Nnode)){
			descL<-smap$edge[which(smap$edge[,1]==node),2][1]
			Lclade<-unique(c(which(smap$edge[,1]==node)[1],which(smap$edge[,2]%in%getDescendants(smap,descL))))
			
			if(!any(sapply(smap$maps[Lclade],function(x)any(names(x)%in%trc)))){ #left descendants DO NOT contain target trait
		
				if(length(Lclade)==1){ 	#if only descendant is tip, drop entire branch
					tbdropped<-c(tbdropped,smap$tip.label[descL])					
				}
		
			}
			
			node = node + 1
			
			}	
		
		smap<-drop.tip.simmap(smap,tbdropped)
		smapL<-smap
					
		tot.len<-length(smap$tip.label)
		node = tot.len + 1
		tbdropped<-c()
		
		while(node <= (tot.len+smap$Nnode)){
			descR<-smap$edge[which(smap$edge[,1]==node),2][2]
			Rclade<-unique(c(which(smap$edge[,1]==node)[2],which(smap$edge[,2]%in%getDescendants(smap,descR))))
			
			if(!any(sapply(smap$maps[Rclade],function(x)any(names(x)%in%trc)))){ #right descendants DO NOT contain target trait
		
				if(length(Rclade)==1){ 	#if only descendant is tip, drop entire branch
					tbdropped<-c(tbdropped,smap$tip.label[descR])
					
				} 
				
			}
			
			node = node +1
			
			}
			
		smap<-drop.tip.simmap(smap,tbdropped)
		
		if(length(smapL$tip.label)==length(smap$tip.label)){
		break
		}
		
		
		}
					
	return(smap)
	}	

