likelihood_subgroup_model<-function(data,phylo,geography.object,model=c("MC","DDexp","DDlin"),par,return.z0=FALSE,maxN=NULL){
 	if(any(grepl("___",phylo$tip.label))){stop("script will not work with '___' in tip labels; remove extra underscores")}
	if(is.null(names(data))){stop("data is unnamed; names should match tip labels in phylogeny")}
	if(!model%in%c("MC","DDexp","DDlin") | length(model)>1){stop("model must be specified as either 'MC', 'DDexp', or 'DDlin'")}
  	if(length(par)!=2){stop("par must contain two values, one for sig2 and another for S, b, or r depending on model")}
	geo.sorted<-.resortGeoObject(phylo,geography.object) 
	
	##create VCV
	if(model=="MC"){
		params<-c(0,par[1],-abs(par[2]))
		mc.ob<-.createModel_MC_geo(phylo,geo.sorted)
        tipdistribution <- getTipDistribution(mc.ob, c(params))            
		V<-tipdistribution$Sigma
	}
	if(model=="DDexp"){
		params<-c(0,par)
		ddexp.ob<-.createModel_DDexp_geo(phylo,geo.sorted)
        tipdistribution <- getTipDistribution(ddexp.ob, c(params))            
		V<-tipdistribution$Sigma
	}
	if(model=="DDlin"){
		params<-c(0,par)
		if(is.null(maxN)){stop("provide maxN value (see help file)")}
		test=exp(par[1])+(par[2]*maxN)
		if(test<=0){return(1000000)}
		ddlin.ob<-.createModel_DDlin_geo(phylo,geo.sorted)
        tipdistribution <- getTipDistribution(ddlin.ob, c(params))            
		V<-tipdistribution$Sigma
	}
	
	##remove rows not in data (i.e., lineages outside of group at present)
	
	itoremove<-which(!colnames(V)%in%names(data))
	if(length(itoremove)>0){V<-V[-itoremove,-itoremove]}
	
	##calculate likelihood
	
	if(any(is.na(V))){
		return(1000000)
		} else{
  		op <- getOption("show.error.messages")
  		options(show.error.messages=FALSE)
		IV=try(solve(V))
  		options(show.error.messages=op)
  		if(class(IV)=="try-error"){
  			IV=pseudoinverse(V) 
  			if(max(IV)==0){return(1000000)}
  		}
 		data<-as.matrix(data[rownames(V)])
		I<-matrix(rep(1,ncol(V)))
		theta<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		if(return.z0==TRUE){return(theta)}
		D<-(I%*%theta)-data[,1]
		T<-t(D)
    	log_final_lik<- -0.5*(T%*%IV%*%D)-0.5*determinant(V)$modulus-0.5*(ncol(V)*log(2*pi))
    	if(is.na(log_final_lik) | is.infinite(log_final_lik)){log_final_lik=-1000000} 
    	return(-log_final_lik) 
    }
}

