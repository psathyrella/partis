likelihood_t_DD_geog<-function(phylo,data,par,geo.object,model=c("DDlin","DDexp"),maxN=NA){
	if(is.null(names(data))){stop("data missing taxa names")}
  	if(length(par)!=2){stop("par must contain two values, one for sig2 and another for the slope")}
	sig2<-exp(par[1])
	rate<-par[2]
	if(model=="DDlin"){
	if(is.na(maxN)){stop("provide maxN value (see help file)")}
	test=sig2+(rate*maxN)
	if(test<=0){return(Inf)}	
	try(V<-.VCV.rescale.DDlin_geog(phylo,sig2,rate,geo.object,check=FALSE))
	} 
	else if(model=="DDexp"){
		try(V<-.VCV.rescale.DDexp_geog(phylo,sig2,rate,geo.object,check=FALSE))
	}
	if(class(V)=="try-error"){return(Inf)}
	if(any(is.na(V))){
		return(Inf)
	} else{
  	op <- getOption("show.error.messages")
  	options(show.error.messages=FALSE)
	IV=try(solve(V))
  	options(show.error.messages=op)
  if(class(IV)=="try-error"){
    IV=corpcor::pseudoinverse(V)
  	if(max(IV)==0){return(Inf)}
  }
  
	data<-as.matrix(data[rownames(V)])
	
	I<-matrix(rep(1,length(phylo$tip.label)))
	theta<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
	D<-(I%*%theta)-data[,1]
	T<-t(D)
	    log_final_lik<- -0.5*(T%*%IV%*%D)-0.5*determinant(V)$modulus-0.5*(length(phylo$tip.label)*log(2*pi))
    
    if(is.na(log_final_lik) | is.infinite(log_final_lik)){log_final_lik=-1000000} 
    return(as.numeric(-log_final_lik)) 
	}
}