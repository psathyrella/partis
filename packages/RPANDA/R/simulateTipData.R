setGeneric(
    name="simulateTipData",
    def=function(object="PhenotypicModel", params="numeric", method="numeric", v="logical"){standardGeneric("simulateTipData")}
)

setMethod(
    f="simulateTipData",
    signature="PhenotypicModel",
    definition=function(object, params, method=3, v=TRUE){
        if(v){
            cat("*** Simulation of tip trait values ***\n")
            beginning <- Sys.time()
        }

        if( method == 1 ){
            if(v){ cat("Computes the tip distribution, and returns a simulated dataset drawn in this distribution.\n") }
            tipdistribution <- getTipDistribution(object, params)
            X <- rmvnorm(1, tipdistribution$mean, tipdistribution$Sigma)
            X <- matrix(data=X, ncol=1)
            rownames(X) <- object@tipLabels

        }else{
            if( method==2 ){
                if(v){ cat("Simulates step-by-step the whole trajectory, plots it, and returns tip data.\n") }
                firstGraph <- TRUE
                # Here you can try "rainbow()", "terrain.colors()", "heat.colors()", "topo.colors()"
                colorChoice <- terrain.colors(length(object@tipLabels))
                xlim <- c(0,object@period[length(object@period)])
                # we get ylim by simulating a few tip distributions and taking the max and min of the results
                findLim <- c(object@initialCondition(params)$mean[1])
                for(k in 1:5){
                    findLim <- c(findLim, simulateTipData(object, params, method=1, v=FALSE) )
                }
                ylim <- c( min(findLim) , max( findLim ) )
                ylim <- c( ylim[1] - (ylim[2]-ylim[1])/10, ylim[2] + (ylim[2]-ylim[1])/10)
            }else{ if(v){ cat("Simulates step-by-step the whole trajectory, but returns only the tip data.\n") } }
            
            initialCondition <- object@initialCondition(params)
            X <- rnorm(length(initialCondition$mean), initialCondition$mean, initialCondition$var)
            dt <- 0.001
            sqrtdt <- sqrt(dt)
        
            for(i in 1:(length(object@period)-1)){
                # If there is a branching event at the beginning of the period
                if(object@numbersPaste[i] != 0){
                    if( object@numbersPaste[i] <= length(X) ){
                        X <- c( X[1:(object@numbersPaste[i]-1)], X[object@numbersCopy[i]], X[object@numbersPaste[i]:length(X)] )
                    }else{
                        X <- c( X, X[object@numbersCopy[i]] )
                    }
                }

                # The evolution of X on the sliced period
                time <- seq( from = object@period[i], to = object@period[i+1], by = dt )

                if( method==2 ){
                    # We remember the value of X at each time step, in the matrix Xt of dimension (nTimeSlices+1)*nLineages
                    nLineages <- length(X)
                    nTimeSlices <- length(time)
                    Xt <- matrix( rep(0, times=(nTimeSlices+1)*nLineages), ncol=nLineages )
                    Xt[1,] <- X
                    for(j in 1:nTimeSlices){
                        aAGammai <- object@aAGamma(i, params)  
                        Xt[j+1,] <- Xt[j,] + (aAGammai$a(time[j]) - aAGammai$A %*% Xt[j,])*dt + sqrtdt* aAGammai$Gamma(time[j]) %*% rnorm(nLineages, 0, 1)
                    }
                    X <- Xt[(nTimeSlices+1),]
                    # Then the trajectories are plotted through time
                    if( firstGraph ){
                        matplot(time, Xt[2:(nTimeSlices+1),], ylab="Trait value", xlab="Time", col=colorChoice[1:nLineages], main="Whole trajectory of trait evolution", cex.main = 1, type="l", lty="solid", xlim=xlim, ylim=ylim)
                        firstGraph <- FALSE
                    }else{
                        matlines(time, Xt[2:(nTimeSlices+1),], lty="solid", col=colorChoice[1:nLineages])
                    }

                }else{
                    # If we don't plot the trajectories, it is faster to forget the previous time step
                    for(t in time){
                        aAGammai <- object@aAGamma(i, params)  
                        X <- X + (aAGammai$a(t) - aAGammai$A %*% X)*dt + sqrtdt* aAGammai$Gamma(t) %*% rnorm(length(X), 0, 1)
                    }
                }
            }

	        X <- matrix(data=X, ncol=1)
            rownames(X) <- object@tipLabelsSimu
        }

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }  

        return(X)
    }
)
