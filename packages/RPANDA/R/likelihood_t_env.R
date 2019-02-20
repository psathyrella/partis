################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################



likelihood_t_env<-function(phylo, data, model=c("EnvExp", "EnvLin"), ...){

# Use mvMORPH for computing the log-likelihood
#require(mvMORPH)

## Parameterization
    par<-list(...)

# Default model
    if(!is.function(model)){
        model<-model[1]
    }

# Number of tips
    tips <- length(phylo$tip.label)
    
# Check if the climatic function is provided
    if(is.null(par[["fun"]])){
        stop("Please provide a time-function")
    }else if(!is.function(par$fun)){
        
        # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
        if(is.null(par[["df"]])){
            par$df <- smooth.spline(par$fun[,1], par$fun[,2])$df
        }
        spline_result <- sm.spline(par$fun[,1],par$fun[,2], df=par$df)
        env_func <- function(t){predict(spline_result,t)}
        
        # if we don't provide a time step in par we take the time steps of the dataset?
        t<-unique(par$fun[,1])
        
        # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
        # the user can choose by specifying it in the "par" list
        
        if(is.null(par[["scale"]])){
            # We build the interpolated smoothing spline function
            par$fun<-splinefun(t,env_func(t))
        }else{
            curve_int<-env_func(t)
            curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
            par$fun<-splinefun(t,curve_scaled)
        }
        
        # Otherwise we assume that the environnemental function is given
    }
    
# Check if the branching time is provided
    if(is.null(par[["times"]])){
        warning("The branching time for the \"phylo\" object was not provided by the user")
        par$times<-branching.times(phylo)
        # root age
        mtot<-max(par$times)
        # Set the root to zero
        par$times<-max(par$times)-par$times
    }else{
# root age
    mtot<-par$mtot
    }
    
# Check if the difference between tip and present day.
    if(is.null(par[["maxdiff"]])){
        maxdiff<-0
    }else{
        maxdiff<-par$maxdiff
    }

# Check if the root value is provided (could be used with an mcmc setting)
    if(is.null(par[["mu"]])){
        par$mu<-NULL
    }

# Check if the tree is in prunning-wise order
    if(is.null(par[["check"]])){
        par$check<-TRUE
    }

# Check if there is polytomies
    if (phylo$Nnode != tips - 1) {
        stop("You can't use this function with polytomies, transform the tree using \"multi2di\" function first")
    }
    
# Check if measurment error is provided
    if(is.null(par[["error"]])){
        is_error<-FALSE
    }else{
        is_error<-TRUE
    }

# Check if measurment error is estimated
    if(is.null(par[["error_param"]])){
      error_param<-FALSE
    }else{
      error_param<-par$error_param
    }

# Control for the integrate function: number of subdivisions
    if(is.null(par[["subdivisions"]])){
        subdivisions<-200L
    }else{
        subdivisions<-par$subdivisions
    }
    
## Transform the tree and return the log-likelihood
    
    # Check the parameters
    if(is.null(par[["param"]]))  {
         stop("Please provide parameters values for \"sig2\" and \"beta\" ")
    }
     
    # Sigma is not provided but analytically computed instead
    phylo <- .CLIMtransform(phylo, param=par$param, mtot=mtot, times=par$times, funEnv=par$fun, model=model, tips=tips, subdivisions=subdivisions, maxdiff=maxdiff)
   
    # Add measurement error
    if(is_error){
        if(error_param){
            phylo$edge.length[par$index_error]<-phylo$edge.length[par$index_error]+par$param[length(par$param)]*par$param[length(par$param)] # estimate the error
        }else{
            phylo$edge.length[par$index_error]<-phylo$edge.length[par$index_error]+par$error^2 # assume the "se" are provided in the error vector
        }
    }
   
 
    # Compute the log-likelihood
    LL <- mvLL(phylo,data,method="pic",param=list(estim=FALSE, check=par$check, mu=par$mu, sigma=1))$logl


if(is.na(LL) | is.infinite(LL)){
return(-1000000)
} # If we use Infinity some optimizer fails; e.g. L-BFGS-B
return(LL)
    
}


##---------------------Functions_used_internally--------------------------------##
## Need to remove it from the export list in the NAMESPACE
## Should I get a general wrapper when neither EnvExp or EnvLin are provided?

## Function to scale the tree to parameters of the climatic model
.CLIMtransform<-function(phy, param, mtot, times, funEnv, model, tips, subdivisions, maxdiff){
    
    # Not yet used (fixed at 0), depends on wether the tree have extant species
    # maxdiff<-0
    res <- phy
    
    if(is.function(model)){
        # user defined function
        f<-function(x){ model((mtot+maxdiff)-x, funEnv, param) }
        
    }else if(model=="EnvExp"){
        # define sigma and beta
        sigma<-exp(param[1])
        beta<-param[2]
        # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
        f<-function(x){sigma*exp(beta*funEnv((mtot+maxdiff)-x))}
        
    }else if(model=="EnvLin"){
        # define sigma and beta
        sigma<-exp(param[1])
        beta<-exp(param[2])
        # Clim-lin function
        f<-function(x){sigma+(beta-sigma)*funEnv((mtot+maxdiff)-x)}
    }
    
    # Transforms the branch-lengths of the tree
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age <- times[phy$edge[i, 1] - tips]
        int <- try(integrate(f,lower=age, upper=age+bl, subdivisions=subdivisions,rel.tol = .Machine$double.eps^0.05), silent = TRUE)
        # Try catch if the integrand is divergent
        if(inherits(int ,'try-error')){
            warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
            integ <- NA_real_
        } else {
            integ <- int$value
        }
        res$edge.length[i] <- integ
    }
    phy<-res
    return(phy)
}

##------------------------------------------------------------------------------##
