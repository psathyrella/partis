CreateGeobyClassObject<-function(phylo,simmap,trim.class,ana.events,clado.events,stratified=FALSE,rnd=5){

trc=trim.class

new.map<-.trimSimmap(simmap,trc)

class.object<-CreateClassObject(new.map)
geo.object<-CreateGeoObject_BioGeoBEARS(full.phylo=phylo,trimmed.phylo=new.map,ana.events,clado.events,stratified=stratified)

geot<-round(geo.object$times,rnd)
clat<-round(class.object$times,rnd)
nodeDist<-sort(unique(c(geot,clat)))
nodeDiff<-diff(nodeDist)
if(any(nodeDiff<= (2*(10^-rnd)))){stop("potential rounding error, two time bins very similar, try changing rnd digits")}
	
nat<-list()

u<-0
y<-0

for(i in 1:length(nodeDiff)){

	if((nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #if timing is the same for both
		u = u+1
		y = y+1
		tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
		tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
		hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
		diag(hold.mat)<-1
		nat[[i]]<-hold.mat
	}
	if((nodeDist[i]%in%geot) && (!nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		y = y+1
		tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
		tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
		hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
		diag(hold.mat)<-1
		nat[[i]]<-hold.mat		
	}
	if((!nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		u = u+1
		tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
		tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
		hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
		diag(hold.mat)<-1
		nat[[i]]<-hold.mat		
	}
	}

	coll.vector<-vector()
	count=1
	geography.matrix<-list()
	for(i in 1:length(nat)){
		int.mat<-nat[[i]] 
		if(count>1 && !is.character(all.equal(nat[[i]],geography.matrix[[count-1]]))){coll.vector<-c(coll.vector,i)} else{	
		geography.matrix[[count]]<-int.mat
		count=count+1
		}
	}
	nodeDist<-nodeDist[which(!nodeDist%in%nodeDist[coll.vector])]
	nodeDist<-nodeDist-min(nodeDist)
	nodeDiff<-diff(nodeDist)	

return(list(map=new.map,geo.object=list(geography.object=geography.matrix,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object
}