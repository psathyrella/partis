setGeneric(
    name="getDataLikelihood",
    def=function(object="PhenotypicModel", data="numeric", params="numeric", v="logical"){standardGeneric("getDataLikelihood")}
)

setMethod(
    f="getDataLikelihood",
    signature="PhenotypicModel",
    definition=function(object, data, params, v=FALSE){
        if(v){
            cat("*** Computing -log( likelihood ) of tip trait data under a given set of parameters ***\n")
            beginning <- Sys.time()
        }

        if(!is.null(rownames(data))){
            data <- data[object@tipLabels,]
        }

        if(object@constraints(params)){
            n <- length(data)
            tipdistribution <- getTipDistribution(object, params)
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

            dataminusXT <- matrix(data - tipdistribution$mean, nrow=1)
            dataminusX <- matrix(data - tipdistribution$mean, ncol=1)

            ProdVectoriel = dataminusXT %*% IV %*% dataminusX

            calcul <-  (ProdVectoriel + determinant(V)$modulus + n*log(2*pi)) /2
			if(is.na(calcul) | is.infinite(calcul)){calcul=-1000000}
        }else{
            calcul <- -Inf
        }

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }
        return(as.numeric(calcul))
    }
)
