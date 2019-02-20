fit_t_comp<-function(phylo,data,model=c("MC","DDexp","DDlin"),pars=NULL,geography.object=NULL){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(is.null(names(data))){stop("data missing taxa names")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
is_tip <- phylo$edge[,2] <= length(phylo$tip.label)
if(sum(diff(phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp cannot be used with ladderized phylogenies')}

if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),0)}
params0<-c(0,pars)

if(is.null(geography.object)){
	if(model=="MC"){
		mc.ob<-.createModel_MC(phylo)
		opt<-fitTipData(mc.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), S = as.numeric(S), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		ddexp.ob<-.createModel_DDexp(phylo)
		opt<-fitTipData(ddexp.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		ddlin.ob<-.createModel_DDlin(phylo)
		opt<-fitTipData(ddlin.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), b = as.numeric(b), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object)){
	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo<-.resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in Marc code
	if(model=="MC"){
		mc.ob<-.createModel_MC_geo(phylo,sgeo)
		opt<-fitTipData(mc.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), S = as.numeric(S), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		ddexp.ob<-.createModel_DDexp_geo(phylo,sgeo)
		opt<-fitTipData(ddexp.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		ddlin.ob<-.createModel_DDlin_geo(phylo,sgeo)
		opt<-fitTipData(ddlin.ob,data,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), b = as.numeric(b), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}
}