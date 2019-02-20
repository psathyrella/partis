################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##                                                                            ##
##                                                                            ##
################################################################################

## Need to be improved. Transform it to C?
sim_t_env<-function(phylo, param, env_data, model, root.value=0, step=0.001, plot=FALSE, ...){
  
  # options
  arg <- list(...)
  if(missing(phylo)) stop("You must provide a phylogenetic tree!!")
  if(inherits(param,"fit_t.env")){
    fun <- param$model
    env_data <- param$env_func
    root.value <- param$root
    param <- param$param
  }else{
    if(missing(param)) stop("You must provide a vector of parameters!!")
    if(missing(model)) stop("You must provide a function for the model!!")
    if(missing(env_data)) stop("You must provide environmental data!!")
    
    # the default function model
    if(!is.function(model)){
      if(model=="EnvExp"){
        # Clim-Exp function
        fun<-function(x, env_data, param){param[1]*exp(param[2]*env_data(x))}
      }else if(model=="EnvLin"){ 
        # Clim-lin function
        fun<-function(x, env_data , param){param[1]+(param[2]-param[1])*env_data(x)}
      }else{
        stop("User defined model must be a function. Otherwise choose the default \"EnvExp\" or \"EnvLin\" models")
      }
    }else{
      fun <- model
    }
    
  }
  
  # Check if the climatic function is provided
if(!is.function(env_data)){
    
    # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
    if(is.null(arg[["df"]])){
      arg$df <- smooth.spline(env_data[,1], env_data[,2])$df
    }
    spline_result <- sm.spline(env_data[,1],env_data[,2], df=arg$df)
    env_func <- function(t){predict(spline_result,t)}
    
    # if we don't provide a time step in par we take the time steps of the dataset?
    t<-unique(env_data[,1])
    
    # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
    # the user can choose by specifying it in the "par" list
    
    if(is.null(arg[["scale"]])){
      # We build the interpolated smoothing spline function
      env_data<-splinefun(t,env_func(t))
    }else{
      curve_int<-env_func(t)
      curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
      env_data<-splinefun(t,curve_scaled)
    }
    
    # Otherwise we assume that the environnemental function is given
  }
  
  ## BM sim temp
  bmSimT<-function(start,n,dt,age,fun_env,param,rootage, env){
    if(n%%dt==0){
      integ<-n%/%dt
      if(integ!=0){
        ageval<-rootage-(cumsum(rep(dt,integ))+age)
        sigmaval<-sapply(1:integ,function(x){ dt*fun_env(ageval[x], env, param)  })
        bmval<-sapply(1:integ,function(x){ rnorm(1,sd=sqrt(sigmaval[x])) })
      }else{
        ageval<-rootage-((n%%dt)+age)
        bmval<-rnorm(1,sd=sqrt(  (n%%dt)*fun_env(ageval, env, param) )   )
      }
      process<-cumsum(c(start,bmval))
    }else{
      integ<-(n%/%dt)
      if(integ!=0){
        ageval<-rootage-(cumsum(c(rep(dt,integ),n%%dt))+age)
        sigmaval<-c(sapply(1:integ,function(x){dt*fun_env(ageval[x], env, param) }),n%%dt*fun_env(ageval[integ+1], env, param) )
        bmval<-sapply(1:integ,function(x){rnorm(1,sd=sqrt(sigmaval[x]))})
      }else{
        ageval<-rootage-((n%%dt)+age)
        bmval<-rnorm(1,sd=sqrt( ((n%%dt))*fun_env(ageval, env, param)  ))
      }
      process<-cumsum(c(start,bmval))
    }
    return(process)
  }
  
  ## discretize time
  bmdiscT<-function(n,dt,age,rootage){
    if(n%%dt==0){
      integ<-n%/%dt
      if(integ!=0){
        sigmaval<-rootage-(cumsum( c(0,rep(dt,integ)) )+age)
      }else{
        sigmaval<-rootage-(c(0,n%%dt)+age)
      }
    }else{
      integ<-(n%/%dt)
      if(integ!=0){
        sigmaval<-rootage-(cumsum(c(rep(dt,integ),n%%dt))+age)#ne pas standardiser par sig puisque juste pour le plot...
      }else{
        sigmaval<-rootage-(c(0,n%%dt)+age)
      }
    }
    return(sigmaval)
  }
  
  # check the phylo order
  tr<-reorder(phylo,"cladewise")
  Ages<-nodeHeights(tr)[,1]
  rootage<-max(nodeHeights(phylo))
  X<-T<-list()
  N<-Ntip(tr)
  for(i in 1:nrow(tr$edge)){
    if(tr$edge[i,1]==(N+1)){
      X[[i]]<-bmSimT(start=root.value, n=tr$edge.length[i], dt=step, age=Ages[i], fun, param, rootage, env_data) # N+1 = element à la racine; root.age =state at the root
    }else{
      parent<-match(tr$edge[i,1],tr$edge[,2]) # recherche le noeud parent de l'element i (peut par ex. etre N+1 si deuxieme noeud sur l'arbre). renvoie la ligne dans edge
      X[[i]]<-bmSimT(start=X[[parent]][length(X[[parent]])], n=tr$edge.length[i], dt=step, age=Ages[i], fun, param, rootage, env_data)
    }
    T[[i]]<-bmdiscT(tr$edge.length[i],dt=step,age=Ages[i],rootage)
  }
  ## Plot function
  plot_process<-function(X,T,...){
    minX<-min(sapply(X,min))
    maxX<-max(sapply(X,max))
    plot(T[[1]][1],X[[1]][1],ylim=c(minX,maxX),xlim=c(0,max(sapply(T,max))),ylab="phenotype",xlab="time")
    for(i in 1:length(X)) lines(T[[i]],X[[i]])
    if(!is.null(arg[["color"]])) cols<-arg$colors
    else cols<-"black"
    return(cols)
  }
  if(plot==TRUE)  plot_process(X,T,arg)
  # retourne les valeurs
  ind<-sapply(1:N,function(x){which(tr$edge[,2]==x)})
  val<-sapply(1:N,function(x){X[[ind[x]]][length(X[[ind[x]]])] } )
  names(val) <- tr$tip.label
# return the trait value
  return(val)
}


