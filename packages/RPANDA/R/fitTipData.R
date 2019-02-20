setGeneric(
    name="fitTipData",
    def=function(object="PhenotypicModel", data="numeric", params0="numeric", GLSstyle="logical", v="logical"){standardGeneric("fitTipData")}
)

setMethod(
    f="fitTipData",
    signature="PhenotypicModel",
    definition=function(object, data, params0=NULL, GLSstyle=FALSE, v=FALSE){
        if(v){
            cat("*** Fit of tip trait data ***\n")
            cat("Finds the maximum likelihood estimators of the parameters, \nreturns the likelihood and the inferred parameters.\n")
            cat("**WARNING** : This function uses the standard R optimizer \"optim\".\nIt may not always converge well.\nPlease double check the convergence by trying\ndistinct parameter sets for the initialisation.\n")
            beginning <- Sys.time()
        }

        n <- length(data)
        if(!is.null(rownames(data))){
            data <- data[object@tipLabels,]
        }

        n <- length(data)

        # If params0 is not given, we use the 'params0' value contained in the model
        if(is.null(params0)){
            params0 <- object@params0
        }
        # In "GLS-style" mode, there is an analytical expression for the first parameter, namely m0
        if(GLSstyle){
            params0 <- params0[2:length(params0)]
        }

        # computing the mean vector and variance matrix for the model, returns -log(likelihood) (a real number)
        toBeOptimized <- function(params){

            if(GLSstyle){paramsPrVerif <- c(0, params)}else{paramsPrVerif <- params}
            if(object@constraints(paramsPrVerif)){

                if(GLSstyle){
                    tipdistribution <- getTipDistribution(object, c(0,params))
                  
		            V<-tipdistribution$Sigma
		            data<-data[rownames(V)]
		  			op <- getOption("show.error.messages")
		  			options(show.error.messages=FALSE)
					IV=try(solve(V))
		  			options(show.error.messages=op)
		  			if(class(IV)=="try-error"){
		    			IV=pseudoinverse(V) 
		  				if(max(IV)==0){return(Inf)}
		  			}
		
                    I<-matrix(rep(1,n))

                    m0 <-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%as.matrix(data)[,1]

                    dataminusXT <- matrix(data - rep(m0, times=n), nrow=1)
                    dataminusX <- matrix(data - rep(m0, times=n), ncol=1)

                    ProdVectoriel = dataminusXT %*% IV %*% dataminusX

                    calcul <-  (ProdVectoriel + determinant(V)$modulus+ n*log(2*pi)) /2
                    params <- c(m0, params)

                }else{
                    calcul <- getDataLikelihood(object, data, params)
                }

            }else{
                calcul <- -Inf
            }

            return(calcul)
        }

        # looking for the argmin of -log(likelihood) (i.e. argmax of likelihood)
        optimisation <- optim(params0, toBeOptimized)
        inferredParams <- optimisation$par
        # In GLS-style, we got all parameters except the first one, 'm0' that we compute through a last call to getTipDistribution
        if(GLSstyle){
            tipdistribution <- getTipDistribution(object, c(0,inferredParams))

		    V<-tipdistribution$Sigma
		  	op <- getOption("show.error.messages")
		  	options(show.error.messages=FALSE)
			IV=try(solve(V))
		  	options(show.error.messages=op)
		  	if(class(IV)=="try-error"){
		    	IV=pseudoinverse(V) 
		  		if(max(IV)==0){return(Inf)}
		  	}
		    data<-data[rownames(V)]
            I<-matrix(rep(1,n))

            m0 <-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%as.matrix(data)[,1]
            inferredParams <- c(m0, inferredParams)
        }
        names(inferredParams) <- object@paramsNames

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(value = optimisation$value, inferredParams = inferredParams, convergence = optimisation$convergence))
    }
)
