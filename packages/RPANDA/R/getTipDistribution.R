setGeneric(
    name="getTipDistribution",
    def=function(object="PhenotypicModel", params="numeric", v="boolean"){standardGeneric("getTipDistribution")}
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicModel",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through ODE resolution ***\n(Method working for any model)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params)
            ai <- aAGammai$a
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma
            n = length(mean)
            # We now need to build the ODE system such that dSigma/dt = -A Sigma - Sigma A + Gamma
            derivativeSigma <- function(t,y,params){
                Sigma = matrix(y,nrow=n)
                dSigma <- -Ai %*% Sigma - t(Sigma) %*% t(Ai) + Gammai(t) %*% t(Gammai(t))
                return(list(dSigma))
            }

            # And we build a second ODE system such that dm/dt = -Ai m + ai
            derivativemean <- function(t,y,params){
                return(list(-Ai %*% y + ai(t)))
            }

            # We update the vectors of means and covariances through their ODE system resolution
            times <- c(object@period[i], object@period[i+1])
            if((object@period[i+1]-object@period[i])> 1e-15 ){
                mean  <- ode(mean, times, derivativemean)[2, 2:(n+1)]
                sigma <- ode(as.vector(Sigma), times, derivativeSigma)[2, 2:(n*n+1)]
            }
            Sigma = matrix(sigma,nrow=n)
        }
        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicACDC",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for an EB process ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        sigma0 <- aAGamma1$Gamma(0)[1,1]
        r <- 2*log(sigma0 / aAGamma1$Gamma(1)[1,1])

        n = length(object@matrixCoalescenceTimes[1,])
        mean <- rep(m0, n)
        Sigma <- v0 + (sigma0**2/(2*r))*(exp(2*r*object@matrixCoalescenceTimes) - 1)

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels
        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicBM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for a BM ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        b <- aAGamma1$a(0)[1]
        sigma <- aAGamma1$Gamma(0)[1,1]

        mean <- m0 + b*diag(object@matrixCoalescenceTimes)
        Sigma <- v0 + sigma**2*object@matrixCoalescenceTimes

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels
        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicDD",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for a DD ***\n")
            beginning <- Sys.time()
        }
        # Warning : parameters must be given in the following order :
        m0 <- params[1]
        v0 <- params[2]
        r <- params[3]
        sigma0 <- params[4]
        N <- length(object@period)

        vectVariances <- rep(0, times=N)
        buffer <- 0
        for(j in 1:(N-1)){
            vectVariances[j] <- buffer
            buffer <- buffer + exp(2*r*object@nLivingLineages[j])*(object@period[j+1] - object@period[j])
        }
        vectVariances[N] <- buffer

        Sigma <- v0 + sigma0^2 * apply(object@matrixCoalescenceJ, c(1,2), function(x){ vectVariances[x] })
        mean <- rep(m0, times=length(Sigma[,1]))

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels
        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicGMM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Analytical computation of tip traits distribution ***\n(Method working for the GMM model only)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        d1 <-params[3]
        d2 <-params[4]
        S <-params[5]
        sigma <-params[6]

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params) 
            delta_t <- object@period[i]-object@period[i+1]
            n1 <- object@n1[i]
            n2 <- object@n2[i]

            exp_SDelta <- exp(S*delta_t)
            alpha <- (exp_SDelta - 1)^2 / 2
            beta <- ( 1 - exp_SDelta^2 ) / 2
            top_part <- matrix( c(rep(alpha/n1, times=n1^2), rep(beta/n2, times=n1*n2)), nrow=n1 )
            bottom_part <- matrix( c(rep(beta/n1, times=n1*n2), rep(alpha/n2, times=n2*n2)), nrow=n2 )
            exp_DeltaA <- diag(exp_SDelta, n1+n2) + rbind(top_part, bottom_part)

            mean <- exp_DeltaA %*% mean + c( rep(-S*delta_t*(d1+d2)/2 + (1- exp_SDelta^2)*(d1-d2)/4, times=n1), rep(-S*delta_t*(d1+d2)/2 + (1- exp_SDelta^2)*(d2-d1)/4,times=n2) )

            calcul_jaune <- sigma^2*(1/n1 + 1/n2)*((1 - exp_SDelta^4)/(16*S) - (1-exp_SDelta^2)/(4*S) - delta_t/4)
            calcul_orange <- sigma^2*(1/n1 + 1/n2)*((exp_SDelta^4 - 1)/(16*S) - delta_t/4)
            top_part <- matrix( c(rep(calcul_jaune, times=n1^2), rep(calcul_orange, times=n1*n2)), nrow=n1 )
            bottom_part <- matrix( c(rep(calcul_orange, times=n1*n2), rep(calcul_jaune, times=n2*n2)), nrow=n2 )
            integrale_AtA <- diag(sigma^2*(1-exp_SDelta^2)/(2*S), n1+n2) + rbind(top_part, bottom_part)

            Sigma <- exp_DeltaA %*% Sigma %*% t(exp_DeltaA) + integrale_AtA
        }

        mean <- matrix(mean,ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicADiag",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through integrated formula ***\n(Method working for models with a constant, A diagonalizable, and Gamma constant)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params) 
            ai <- aAGammai$a(0)
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma(0)
            n = length(Ai[,1])

            # The mean and covariances evolve through the period with the following formula
            e <- eigen(Ai, symmetric=TRUE)
            Q <- e$vectors
            Qt <- t(Q)

            delta_t <- object@period[i+1]-object@period[i]
            R1 <- diag( exp( -delta_t * e$values ) )
            r2 <- rep(0, n)
            r3 <- rep(0, n)
            for(k in 1:n){ 
                if(abs(delta_t*e$values[k]) > 1e-3 ){
                    r2[k] <- (1 - exp( -delta_t*e$values[k] ))/e$values[k]
                    r3[k] <- (1 - exp( -delta_t*2*e$values[k] ))/(2*e$values[k])
                }else{ # DL en 0 ordre 3 
                    r2[k] <- delta_t + 1/2*e$values[k]*delta_t^2 + 1/6*e$values[k]^2*delta_t^3
                    r3[k] <- delta_t +     e$values[k]*delta_t^2 + 2/3*e$values[k]^2*delta_t^3
                }      
            }
            R2 <- diag( r2 )
            R3 <- diag( r3 )

            mean <- Q %*% R1 %*% Qt %*% mean + Q %*% R2 %*% Qt %*% ai
            Sigma <- Q %*% R1 %*% Qt %*% Sigma %*% Q %*% R1 %*% Qt + Gammai %*% t(Gammai) %*% Q %*% R3 %*% Qt            
        }
           
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicPM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Analytical computation of tip traits distribution ***\n(Method working for models with a constant, A = aI + bU, and Gamma constant)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Works for models with or without the "OU" part
        if( object@aAGamma(1, params)$OU == TRUE ){
            theta <-params[3]
            psi   <-params[4]
            S     <-params[5]
            sigma <-params[6]
        }else{
            theta <-0
            psi   <-0
            S     <-params[3]
            sigma <-params[4]
        }

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params) 
            vectorU <- aAGammai$u
            delta_t <- object@period[i+1]-object@period[i]

            n <- sum(vectorU)
            d <- exp(-(psi+S)*delta_t*vectorU)      # d is a vector
            SigmaU <- as.vector(Sigma%*%vectorU)    # SigmaU is a vector. Note that SigmaU = USigma because A is symmetric.
            f <- (exp(-psi*delta_t) - exp(-(psi+S)*delta_t))/n      # f is a scalar
            Sigma <- outer(d,d)*Sigma + outer(f*d*SigmaU, vectorU) + outer(f*vectorU, SigmaU*d) + outer(f^2 * c(vectorU%*%SigmaU) * vectorU,vectorU) 
            mean <- f*( vectorU%*%mean )*vectorU + d*mean

            # The exponential or the Taylor expansion depending on the parameter values
            if( sigma != 0 ){
                if (abs((psi+S)*delta_t) > 1e-3){
                    integral1 <- (1-exp(-2*(psi+S)*delta_t))/(2*(psi+S))
                }else{
                    integral1 <- delta_t - (psi+S)*delta_t^2
                } 
                if (abs(psi*delta_t) > 1e-3) {
                    integral2 <- (1-exp(-2*psi*delta_t))/(2*psi)
                }else{
                    integral2 <-  delta_t - psi*delta_t^2
                }
                Sigma <- Sigma + diag(sigma^2 * integral1 *vectorU) + outer(sigma^2 * (integral2 - integral1)/n * vectorU, vectorU)
            }

            # No problem for the integral associated with a_i and the evolution of the mean
            mean <- mean + theta * (1-exp(-psi*delta_t)) * vectorU
        }

        mean <- matrix(mean,ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicOU",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for an OU process ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        psi <- aAGamma1$A[1,1]
        theta <- aAGamma1$a(0)[1]/psi
        sigma <- aAGamma1$Gamma(0)[1,1]

        mean <- theta + (m0 - theta)*exp(-psi*diag(object@matrixCoalescenceTimes))
        n = length(mean)
        Sigma <- diag(0, n)
        for(k in 1:n){
            for(l in k:n){
                valeur <- exp(-psi*(object@matrixCoalescenceTimes[k,k] + object@matrixCoalescenceTimes[l,l] - 2*object@matrixCoalescenceTimes[k,l])) * ( sigma**2/(2*psi) + exp(-2*psi*object@matrixCoalescenceTimes[k,l])*(v0 - sigma**2/(2*psi) ) )
                Sigma[k,l] <- valeur
                Sigma[l,k] <- valeur
            }
        }

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels
        return(list(mean = mean, Sigma = Sigma))
    }
)




updateBranchingMatrixSigma <- function(Sigma, copy, paste){
    # copy of a branching lineage in the matrix of covariances 'Sigma'
    n = length(Sigma[1,])
    newSigma <- diag(0, n+1)
    
    newSigma[1:(paste-1),1:(paste-1)] <- Sigma[1:(paste-1),1:(paste-1)]
    newSigma[paste,1:(paste-1)] <- Sigma[copy,1:(paste-1)]
    newSigma[1:(paste-1),paste] <- Sigma[1:(paste-1),copy]
    newSigma[paste,paste] <- Sigma[copy,copy]

    if(paste < n+1){
        newSigma[(paste+1):(n+1),1:(paste-1)] <- Sigma[paste:n,1:(paste-1)]
        newSigma[(paste+1):(n+1),paste] <- Sigma[paste:n,copy]
        newSigma[1:(paste-1),(paste+1):(n+1)] <- Sigma[1:(paste-1),paste:n]
        newSigma[paste,(paste+1):(n+1)] <- Sigma[copy,paste:n]
        newSigma[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma[paste:n, paste:n]
    }

    return(newSigma)
}

