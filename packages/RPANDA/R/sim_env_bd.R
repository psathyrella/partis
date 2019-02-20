sim_env_bd <- function (env_data, f.lamb, f.mu, lamb_par, mu_par, df=NULL, time.stop = 0, return.all.extinct=TRUE, prune.extinct=TRUE)

{

birthdeath.tree.timevar_simp <- function (f.lamb, f.mu, lamb_par, mu_par, time.stop = 0, return.all.extinct=TRUE, prune.extinct=TRUE)

{
    if (time.stop == 0) 
    stop("Must have stopping time")
    
    while (1) {
            	nblineages<-c(1)
            	times<-c(0)
            	b<-f.lamb(0,lamb_par)
            	d<-f.mu(0,mu_par)
            	dt <- rexp(1,(b + d))
            	#print(dt)
            	t<-dt
    			
    			if (t >= time.stop) 
    			{
    					t <- time.stop
    					alive<-1
    					times<-c(times,t)
    					nblineages<-c(nblineages,1)
    					break}
            	
            	r <- runif(1)
            	
            	if (r>b/(b + d))
            	{
            		#print("die")
            		times<-c(times,dt)
            		nblineages<-c(nblineages,0)
            		alive<-rep(FALSE,1)}
            	
            	else
            	{
    				edge <- rbind(c(1, 2), c(1, 3))
    				edge.length <- rep(NA, 2)
    				stem.depth <- rep(t, 2)
    				alive <- rep(TRUE, 2)
    				times<-c(times,dt)
    				nblineages<-c(nblineages,sum(alive))
    				next.node <- 4
    				repeat 
    					{
    					if (sum(alive) == 0)	 
    					break
    					b<-f.lamb(t,lamb_par)
        	    		d<-f.mu(t,mu_par)
    					dt <- rexp(1, sum(alive) * (b + d))
    					t <- t + dt
    					if (t >= time.stop) 
    						{
    						t <- time.stop
    						times<-c(times,t)
    						nblineages<-c(nblineages,sum(alive))
    						break}
    				
    				r <- runif(1)
    				if (r <= b/(b + d)) 
    					{
    					#print("speciation")
    					random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    					e <- matrix(edge[alive, ], ncol = 2)
    					parent <- e[random_lineage, 2]
    					alive[alive][random_lineage] <- FALSE
    					edge <- rbind(edge, c(parent, next.node), c(parent,next.node + 1))
    					next.node <- next.node + 2
    					alive <- c(alive, TRUE, TRUE)
    					stem.depth <- c(stem.depth, t, t)
    					x <- which(edge[, 2] == parent)
    					edge.length[x] <- t - stem.depth[x]
    					edge.length <- c(edge.length, NA, NA)
    					times<-c(times,t)
    					nblineages<-c(nblineages,sum(alive))}
    				
    				else 
    					{
    					#print("extinction")
    					random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    					edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
    					alive[alive][random_lineage] <- FALSE
    					times<-c(times,t)
    					nblineages<-c(nblineages,sum(alive))}  					
    					}
    				 }
    			
    			if (return.all.extinct == TRUE | sum(alive) > 0) 
    				{
    								#print("return.tree")
    								break}
    			}
    							
    							
     							
	if (sum(alive)==0) {obj<-NULL}
 	else if (sum(alive)==1) {obj<-1}
 	else
 	{
 	edge.length[alive] <- t - stem.depth[alive]
 	n <- -1
 	for (i in 1:max(edge)) {
 	if (any(edge[, 1] == i)) {
 	edge[which(edge[, 1] == i), 1] <- n
 	edge[which(edge[, 2] == i), 2] <- n
 	n <- n - 1}}
    edge[edge > 0] <- 1:sum(edge > 0)
    tip.label <- 1:sum(edge > 0)
    mode(edge) <- "character"
    mode(tip.label) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)
    class(obj) <- "phylo"
    obj <- old2new.phylo(obj)
     if (prune.extinct)
    {obj<-drop.extinct(obj)}
    }
    return(list("tree"=obj,"times"=times,"nblineages"=nblineages))
}

if (is.null(df))
  {
    df <- smooth.spline(x=env_data[,1], env_data[,2])$df
  }
  spline_result <- smooth.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  
  lower_bound_control <- 0.10
  upper_bound_control <- 0.10
  lower_bound <- min(env_data[,1])
  upper_bound <- max(env_data[,1])
  # Tabulation of the function from lower_bound -10%, upper_bound + 10%
  time_tabulated <- seq(from=lower_bound*(1.0-lower_bound_control),
                        to=upper_bound*(1.0+upper_bound_control),
                        length.out=1+1e6)
  env_tabulated <- env_func(time_tabulated)
  # Tabulated function
  env_func_tab <- function(t)
  {
    b <- upper_bound * (1.0 + upper_bound_control)
    a <- lower_bound * (1.0 - lower_bound_control)
    # number of intervals
    n <- length(env_tabulated) - 1
    index <- 1 + as.integer( (t - a) * n / (b - a))
    return(env_tabulated[index])
  }
  f.lamb.env <- function(t,y){ f.lamb(t, env_func(time.stop-t)$y, y)}
  f.mu.env <- function(t,y){ f.mu(t, env_func(time.stop-t)$y, y)}
 
 res <- birthdeath.tree.timevar_simp(f.lamb.env, f.mu.env, lamb_par, mu_par, time.stop, return.all.extinct, prune.extinct) 
 
 return(res)
 
 }
