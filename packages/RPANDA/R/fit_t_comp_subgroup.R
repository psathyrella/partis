fit_t_comp_subgroup<-function(full.phylo,ana.events,clado.events,stratified=FALSE,map,data,trim.class,model=c("MC","DDexp","DDlin"),par=NULL,method="Nelder-Mead",bounds=NULL){

	if(is.null(names(data))){stop("data missing taxa names")}
	if(!is.null(dim(data))){stop("data needs to be a single trait")}
	is_tip <- full.phylo$edge[,2] <= length(full.phylo$tip.label)
	if(sum(diff(full.phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_subgroup cannot be used with ladderized full.phylogenies')}
	
	if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }
    
	GeoByClassObject<-CreateGeobyClassObject(full.phylo,map,trim.class,ana.events,clado.events,stratified=stratified)

	phylo<-GeoByClassObject$map
	
	root.trimmed.phylo<-max(nodeHeights(phylo))
	root.data<-max(nodeHeights(drop.tip.simmap(phylo,phylo$tip.label[which(!phylo$tip.label%in%names(data))])))
	
	if(round(root.trimmed.phylo,5)!=round(root.data,5)){stop("error where root of trimmed simmap and root of target clade don't match")}
	#if this above error is triggered, need to use scripts written by JPD using JC mvMORPH approach in summer 2018
	
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	geo.object<-GeoByClassObject$geo.object
	#geo.sorted<-.resortGeoObject(phylo,geo.object) 

	if(length(geo.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(length(data)>length(phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
	if(!all(names(data)%in%phylo$tip.label)){stop("error: some tips missing from pruned simmap")}

	if(is.null(par)){par<-c(log(sqrt(var(data)/max(nodeHeights(extract.clade(phylo,getMRCA(phylo,names(data))))))),0)}
	
	if(model=="MC"){
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="MC",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		S = -abs(opt$par[2])
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="MC",par=opt$par,return.z0=TRUE)
		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		r = opt$par[2]
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDexp",par=opt$par,return.z0=TRUE)
		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		geography.matrix<-geo.object$geography.object
		maxN<-max(vapply(geography.matrix,function(x)max(rowSums(x)),1))
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,maxN=maxN,data=data,model="DDlin",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		b = opt$par[2]
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDlin",par=opt$par,return.z0=TRUE,maxN=maxN)
		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}
