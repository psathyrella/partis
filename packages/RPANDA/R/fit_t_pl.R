################################################################################
##                                                                            ##
##  RPANDA : Leave-One-Out Cross-Validation for High-dimentional              ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################


# Y        = traits matrix (columns: variables, rows: species)
# tree     = phylogenetic tree (an object of class 'phylo')
# model    = either "BM", "OU", "EB", or "lambda"; the model of traits evolution
# method   = "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator. "RidgeAltapprox" and "LASSOapprox" are fast approximations of the LOOCV for the Ridge quadratic and LASSO penalties
# targM    = "null", "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt", and "RidgeAltapprox" methods.
# REML     = TRUE (default) or FALSE. The likelihood method used to estimate the parameters. REML must be preferred with small sample size in order to obtain unbiased estimates.
# up       = upper bound for model parameter search
# low      = lower bound for model parameter search
# tol      = lower bound tolerance for the regularizer (tuning parameter)
# starting = starting values for the parameter search. Must be a vector likes: c(model, regularization)


#require(mvMORPH)    # >= 1.1.0
#require(glassoFast) # https://github.com/JClavel/glassoFast

fit_t_pl <- function(Y, tree, model=c("BM","OU","EB","lambda"), method=c("RidgeAlt","RidgeArch","RidgeAltapprox","LASSO","LASSOapprox"), targM=c("null","Variance","unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=NULL, starting=NULL, SE=NULL, scale.height=TRUE, ...){
    
    # Preliminary checks
    if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
    if(nrow(Y) != Ntip(tree)) stop("Length of phenotypic and phylogenetic data do not match")
    if (all(rownames(Y) %in% tree$tip.label)){
        Y <- Y[tree$tip.label,]
    }else{
        warning("rownames in Y are missing. It is assumed that they are in the same order as in the tree.")
    }
    
    # Supp. param
    parLoo <- list(...)
    if(is.null(parLoo[["echo"]])){ echo <- TRUE }else{ echo <- parLoo$echo}
    
    # Select the model
    model <- match.arg(model)[1]
    
    # Select the method
    method <- match.arg(method)[1]
    
    # Select the method
    targM <- match.arg(targM)[1]
    
    # Bounds for models
    if(is.null(up)){
        switch(model,
        "EB"={up <- 0},
        "OU"={up <- 10},
        "lambda"={up <- 1.1})
    }
    
    if(is.null(low)){
        switch(model,
        "EB"={low <- -10},
        "OU"={low <- 0},
        "lambda"={low <- 1e-5})
    }
    
    # Parameters
    if(ncol(Y)==1) stop("Only works with multivariate datasets")
    n <- nO <- nrow(Y)
    if(model=="OU" & !is.ultrametric(tree)) nC <- n else nC <- n-1
    p <- ncol(Y)
    if(REML==TRUE) n <- n-1
    
    # Scale the tree to unit length
    if(scale.height==TRUE) tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
    
    # Identity matrix
    I <- diag(p)
    
    # Default design for contrasts
    X <- matrix(1, ncol=1, nrow=nO)
    
    # Default penalty is Ridge "null"
    if(method=="RidgeArch" & targM=="null"){
        warning("The \"null\" target cannot be used with the \"RidgeArch\" method. The \"unitVariance\" target is used instead.")
        targM <- "unitVariance"
    }else{
        target <- matrix(0,p,p)
    }
    
    # Warning for rotation invariance
    if(method=="LASSO" | method=="LASSOapprox" | targM=="Variance") warning("The LASSO penalty and \"Variance\" target should not be used on rotated data (see Clavel et al. 2018 for further details on applicability)")
    
    # Identifying tips values
    tipsIndices <- which(tree$edge[, 2] <= Ntip(tree))
    if(!is.null(SE) && SE!=TRUE) SE <- NULL
    
    
    # Default tolerance for the parameter search
    if(is.null(tol)){
        if(method=="RidgeArch"){
            tol = 1e-8
        }else{
            tol = 0
        }
    }
    
    ## ---- Set parameters and bounds
    switch(model,
    "lambda"={
        
        if(method!="RidgeArch"){
            upperBound <- c(log(10e6), up)
            lowerBound <- c(log(tol), low)
        }else{
            upperBound <- c(1, up)
            lowerBound <- c(tol, low)
        }
        
        idx1 <- 1; idx2 <- 2; idx3 <- 3
    },
    "OU"={
        
        if(method!="RidgeArch"){
            upperBound <- c(log(10e6), log(up))
            lowerBound <- c(log(tol), log(low))
        }else{
            upperBound <- c(1, log(up))
            lowerBound <- c(tol, log(low))
        }
        
        idx1 <- 1; idx2 <- 2; idx3 <- 3
    },
    "EB"={
        
        if(method!="RidgeArch"){
            upperBound <- c(log(10e6), up)
            lowerBound <- c(log(tol), low)
        }else{
            upperBound <- c(1, up)
            lowerBound <- c(tol, low)
        }
        
        idx1 <- 1; idx2 <- 2; idx3 <- 3
    },
    "BM"={

        if(method!="RidgeArch"){
            upperBound <- log(10e6)
            lowerBound <- log(tol)
        }else{
            upperBound <- 1
            lowerBound <- tol
        }
        
        idx1 <- idx2 <- 1
        idx3 <- 2
    })
    
    ## ------ Leave-One-Out Cross-validation / Penalty

    loocv <- function(par){
            
            # parameters
            mod_par = par[idx2]
            
            if(!is.null(SE)) error = par[idx3]*par[idx3] else error = NULL
            # Transform the tree
            corrstruct <- .transformTree(tree, mod_par, model=model, mserr=error, Y=Y, X=X, REML=REML)
            if(any(!is.finite(corrstruct$Y))) return(1e6)
            B <- pseudoinverse(corrstruct$X)%*%corrstruct$Y
            residuals <- corrstruct$Y - corrstruct$X%*%B
            
            # Covariance
            Sk <- crossprod(residuals)/n
     
           switch(method,
           "LASSOapprox"={
               # Regularization parameter
               alpha = exp(par[idx1])
               
               # LASSO penalty
               LASSO <- glassoFast(Sk, alpha, maxIt=5000)
               G <- LASSO$w
               Gi <- LASSO$wi
               
               # LOO cross-validated log-likelihood
               LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
                   # Indices matrice
                   Ind <- ifelse(Iridge==0,0,1);
                   
                   Tk <- sapply(1:nC, function(i){
                       Sk <- tcrossprod(Y[i,]); ;
                       A <-.vec((Sridge - Sk)*Ind);
                       BC <- Iridge%*%((Semp - Sk)*Ind)%*%Iridge;
                       sum(A*BC)
                   })
                   
                   bias <- (1/(2*n*(n-1))) * sum(Tk)
                   
                   return(bias)
               }
               
               klbias <- LOObias(Sk, Gi, G, residuals, alpha)
               
               ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*corrstruct$det + n*determinant(G)$modulus + n*sum(Gi*Sk))) + klbias
           },
           "LASSO"={
               # Regularization parameter
               alpha = exp(par[idx1])
               
               # log-lik
               llik <- sapply(1:nC, function(x){
                   Sk <- crossprod(residuals[-x,])/(n-1);
                   LASSO <- glassoFast(Sk, alpha, maxIt=500);
                   G <- LASSO$w;
                   Gi <- LASSO$wi;
                   Swk <- tcrossprod(residuals[x,]);
                   rk <- sum(Swk*Gi);
                   determinant(G)$modulus + rk
               })
               
               
               ll <- 0.5 * (n*p*log(2*pi) + p*corrstruct$det + n*mean(llik))
           },
           "RidgeArch"={
               # Regularization parameter
               alpha = par[idx1]
               
               # Target
               if(targM=="Variance"){
                   target <- diag(diag(Sk))
               }else if(targM=="unitVariance"){
                   target <- I*mean(diag(Sk))
               }
               
               # Hoffbeck & Landgrebe 1996 parameterization
               beta <- (1 - alpha)/(n - 1)
               G <- n*beta*Sk + alpha * target
               Gi <- try(chol(G), silent = TRUE)
               if(inherits(Gi, 'try-error')) return(1e6)
               
               # log-lik
               llik <- sapply(1:nC, function(x){
                   # log-lik form of Hoffbeck & Landgrebe 1996
                   rk <- sum(backsolve(Gi, residuals[x,], transpose = TRUE)^2)
                   log(1 - beta*rk) + (rk/(1 - beta*rk))
               })
               
               ll <- 0.5 * (n*p*log(2*pi) + p*corrstruct$det + n*sum(2*log(diag(Gi))) + n*mean(llik))
           },
           "RidgeAlt"={
               # Regularization parameter
               alpha = exp(par[idx1])
               
               # Choose an updated target matrix if it's not the usual Ridge
               if(targM=="Variance"){
                   target <- diag(1/diag(Sk))
               }else if(targM=="unitVariance"){
                   target <- I*(1/mean(diag(Sk)))
               }
               
               # log-lik
               llik <- try(sapply(1:nC, function(x){
                   Sk <- crossprod(residuals[-x,])/(n-1)
                   pen <- .makePenalty(Sk,alpha,target,targM)
                   Gi <- pen$P
                   detG <- sum(log(pen$ev))
                   Swk <- tcrossprod(residuals[x,])
                   rk <- sum(Swk*Gi)
                   detG + rk
               }), silent = TRUE)
               
               if(inherits(llik, 'try-error')) return(1e6)
               # det of the phylo matrix
               ll <- 0.5 * (n*p*log(2*pi) + p*corrstruct$det + n*mean(llik))
           },
           "RidgeAltapprox"={
               # Regularization parameter
               alpha = exp(par[idx1])
               
               # Choose an updated target matrix if it's not the usual Ridge
               if(targM=="Variance"){
                   target <- diag(1/diag(Sk))
               }else if(targM=="unitVariance"){
                   target <- I*(1/mean(diag(Sk)))
               }
               
               # Ridge penalty
               G <- .makePenalty(Sk,alpha,target,targM)$S
               if (any(!is.finite(G))) return(1e6)
               eig <- eigen(G, symmetric=TRUE) # we can return the vectors and inverse from .makePenalty directly
               V <- eig$vectors
               d <- eig$values
               Gi <- V%*%diag(1/d)%*%t(V)
               H <- (1/(kronecker(d,d)+alpha))
               
               
               # LOO cross-validated log-likelihood
               LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
                   Tk <- sapply(1:nC, function(i){
                       Sk <- tcrossprod(Y[i,]);
                       VSV <- .vec(t(V)%*%((Sridge - (Sk - lambda*target) - lambda*Iridge))%*%(V));
                       sum(VSV * H*.vec(t(V)%*%(-I%*%(Sk - Semp))%*%(V)))
                   })
                   
                   bias <- (1/(2*n*(n-1))) * sum(Tk)
                   
                   return(bias)
               }
               
               klbias <- LOObias(Sk, Gi, G, residuals, alpha)
               
               ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*corrstruct$det + n*sum(log(d)) + n*sum(Gi*Sk))) + klbias
           })
            if (!is.finite(ll)) return(1e6)
            return(ll)
    }
    
    
    # Starting values over range of parameters for quadratic ridge and LASSO (because the tuning value is between 0->Inf)
    # computationally intensive but maybe better to ensure good starting values
    
    if(is.null(starting)){
        
        # Here we can use the forking to spread the calculus over the grid on several cores
        # we can use a randomized search to speed up the computations of very complex models
        
        if(!is.null(SE)){
            # various errors
            guess <- c(0.001,0.01,0.1,1,10)
            error_guess = sqrt(guess)
            lowerBound = c(lowerBound,0)
            upperBound = c(upperBound,Inf)
        }
        
       if(echo==TRUE) message("Initialization via grid search. Please wait...")
        if(method=="RidgeArch"){
            range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
        }else{
            range_val <- log(c(1e-12, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
        }
        switch(model,
        "lambda"={
            mod_val <- c(0.2,0.5,0.8)
            if(!is.null(SE)){
                brute_force <- expand.grid(range_val,mod_val,error_guess)
            }else{
                brute_force <- expand.grid(range_val,mod_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[1]
        },
        "OU"={
            mod_val <- log(log(2)/(max(node.depth.edgelength(tree))/c(0.1,0.5,1.5,3,8)))
            if(!is.null(SE)){
                brute_force <- expand.grid(range_val,mod_val,error_guess)
            }else{
                brute_force <- expand.grid(range_val,mod_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[1]
        },
        "EB"={
            mod_val <- -log(2)/(max(node.depth.edgelength(tree))/c(0.1,0.5,1.5,3,8))
            if(!is.null(SE)){
                brute_force <- expand.grid(range_val,mod_val,error_guess)
            }else{
                brute_force <- expand.grid(range_val,mod_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[1]
        },
        "BM"={
            mod_val = NULL
            if(!is.null(SE)){
                brute_force <- expand.grid(range_val,error_guess)
            }else{
                brute_force <- expand.grid(range_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[1]
        })
        if(echo==TRUE){
            if(method=="RidgeArch"){
                cat("Best starting for the tuning: ",as.numeric(tuning))
            }else{
                cat("Best starting for the tuning: ",as.numeric(exp(tuning)))
            }
        }
    }else{
        start <- .starting_val(starting, model, method, SE)
    }
    # Initial guesses found we start the optimization
    if(echo==TRUE) message("Start optimization. Please wait...")
    
    # Optimization of the cross-validated likelihood
    estimModel <- optim(start, fn = loocv, method="L-BFGS-B", upper=upperBound, lower=lowerBound)
    
    # Compute the scaled tree
    if(!is.null(SE)) SE = (estimModel$par[idx3])*(estimModel$par[idx3])
    corrstruct <- .transformTree(tree, estimModel$par[idx2], model=model, mserr=SE, Y=Y, X=X, REML=REML)
    B <- pseudoinverse(corrstruct$X)%*%corrstruct$Y
    residuals <- corrstruct$Y - corrstruct$X%*%B
    
    # Estimated value for the model parameter
    switch(model,
    "BM"={ model.par <- 0 },
    "OU"={ model.par <- exp(estimModel$par[idx2])},
     model.par <- estimModel$par[idx2]
    )
    
    # Estimated value for the regularization parameter
    if(method=="RidgeArch"){
       gamma <- estimModel$par[idx1]
    }else{
       gamma <- exp(estimModel$par[idx1])
    }
    
    # Compute R
    Snew <- crossprod(residuals)/n
    matMeth <- method
    if(method == "LASSOapprox"){
        matMeth <- "LASSO"
    }else if(method == "RidgeAltapprox"){
        matMeth <- "RidgeAlt"
    }
    regularizedEstimates <- .covPenalized(S=Snew, method=matMeth, targM=targM, tuning=gamma)
    variables <- list(Y=Y, tree=tree)
    # End
    if(echo==TRUE) message("Done in ", estimModel$count[1]," iterations.")
    
    # return the results
    results <- list(loocv=estimModel$value, model.par=model.par, gamma=gamma, corrstruct=corrstruct, model=model, method=method, p=p, n=nO, targM=targM, R=regularizedEstimates, REML=REML, SE=SE, variables=variables)
    class(results) <- "fit_pl.rpanda"
    return(results)
}

## -------------- Miscellaneous functions

# Alternative penalty of van Wieringen & Peeters 2016 - Computational Statistics and Data Analysis
# see also Witten & Tibshirani 2009
.makePenalty <- function(S, lambda, target, targM){
    
    switch(targM,
    "Variance"={
        D <- (S - lambda * target)
        D2 <- D %*% D
        sqrtM <- .sqM(D2/4 + lambda * diag(nrow(S)))
        Alt <- D/2 + sqrtM
        AltInv <- (1/lambda)*(Alt - D)
        evalues <- eigen(Alt, symmetric=TRUE, only.values = TRUE)$values
    },
    "unitVariance"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values - lambda*target[1]
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- diag(evalues)
        D2 <- diag(1/evalues) # Inverse
        Alt <- Q %*% D1 %*% t(Q)
        AltInv <- Q %*% D2 %*% t(Q)
    },
    "null"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- diag(evalues)
        D2 <- diag(1/evalues)
        Alt <- Q %*% D1 %*% t(Q)
        AltInv <- Q %*% D2 %*% t(Q)
    }
    )
    pen <- list(S=Alt, P=AltInv, ev=evalues)
    return(pen)
}

# Matrix square root
.sqM <- function(x){
    if(!all(is.finite(x))) return(Inf)
    eig <- eigen(x, symmetric = TRUE)
    if(any(eig$values<0)) eig$values <- abs(eig$values) # to avoid warning message?
    sqrtM <- tcrossprod(eig$vectors %*% diag(sqrt(eig$values)), eig$vectors)
    return(sqrtM)
}

# Compute the Regularized covariance and it's inverse
.covPenalized <- function(S, method, targM="null", tuning=0){
    
    # dim of S
    p = ncol(S)
    
    # init the target matrix
    if(method=="RidgeAlt"){
        switch(targM,
        "null"={Target <- matrix(0,p,p)},
        "Variance"={Target <- diag(1/diag(S))},
        "unitVariance"={Target <- diag(1/mean(diag(S)),p)})
    }else if(method=="RidgeArch"){
        switch(targM,
        "null"={Target <- matrix(0,p,p)},
        "Variance"={Target <- diag(diag(S))},
        "unitVariance"={Target <- diag(mean(diag(S)),p)})
    }
    
    # Construct the penalty term
    switch(method,
    "RidgeAlt"={
        pen <- .makePenalty(S,tuning,Target,targM)
        P <- pen$S
        Pi <- pen$P
    },
    "RidgeArch"={
        P <- (1-tuning)*S + tuning*Target
        eig <- eigen(P)
        V <- eig$vectors
        d <- eig$values
        Pi <- V%*%diag(1/d)%*%t(V)
    },
    "LASSO"={
        LASSO <- glassoFast(S,tuning)
        P <- LASSO$w
        Pi <- LASSO$wi
    },
    "ML"={
        P <- S
        eig <- eigen(P)
        V <- eig$vectors
        d <- eig$values
        Pi <- V%*%diag(1/d)%*%t(V)
    }
    )
    estimate <- list(R=P, Rinv=Pi)
    return(estimate)
}


# Vec operator
.vec <- function(x) as.numeric(x)

# Trace operator
.tr <- function(x) sum(diag(x))

# Transform starting values
.starting_val <- function(starting, model, method, SE=NULL){
    switch(model,
    "OU"={
        if(method=="RidgeArch"){
            start <- c(starting[1], log(starting[2]))
        }else{
            start <- c(log(starting[1]), log(starting[2]))
        }
    },
    "EB"={
        if(method=="RidgeArch"){
            start <- c(starting[1], starting[2])
        }else{
            start <- c(log(starting[1]), starting[2])
        }
    },
    "BM"={
        if(method=="RidgeArch"){
            start <- starting[1]
        }else{
            start <- log(starting[1])
        }
    },
    "lambda"={
        if(method=="RidgeArch"){
            start <- c(starting[1], starting[2])
        }else{
            start <- c(log(starting[1]), starting[2])
        }
    })
    
    if(!is.null(SE) & model!="BM") start <- c(start, sqrt(starting[3]))
    if(!is.null(SE) & model=="BM") start <- c(start, sqrt(starting[2]))
    return(start)
}


# ------------------------------------------------------------------------- #
# .transformTree                                                            #
# options: phy, param, model, mserr=NULL, Y=NULL, X=NULL, REML=TRUE         #
#                                                                           #
# ------------------------------------------------------------------------- #

.transformTree <- function(phy, param, model=c("EB", "BM", "lambda", "OU"), mserr=NULL, Y=NULL, X=NULL, REML=TRUE){
    
    # pre-compute and checks
    if(attr(phy,"order")!="cladewise") phy <- reorder.phylo(phy, "cladewise")
    n <- Ntip(phy)
    parent <- phy$edge[,1]
    descendent <- phy$edge[,2]
    extern <- (descendent <= n)
    N <- 2*n-2
    diagWeight <- NULL
    
    # Model
    switch(model,
    "OU"={
        param = exp(param)
      if(param!=0){
        D = numeric(n)
        
        # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
        if(!is.ultrametric(phy)){
            dis = node.depth.edgelength(phy) # has all nodes
            D = max(dis[1:n]) - dis[1:n]
            D = D - mean(D)
            phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
        }
        
        # Branching times (now the tree is ultrametric)
        times <- branching.times(phy)
        Tmax <- max(times)
        # compute the branch lengths
        distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
        d1 = distRoot[parent-n]
        d2 = numeric(N)
        d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
        d2[!extern] = distRoot[descendent[!extern]-n]
        
        # weights for a 3 points structured matrix
        diagWeight = exp(param*D)
        phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
        names(diagWeight) = phy$tip.label
        
        # transform the variables
        w <- 1/diagWeight
        Y <- matrix(w*Y, nrow=n)
        X <- matrix(w*X, nrow=n)
        
        # Adjust errors
        if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
      }
    },
    "EB"={
        if (param!=0){
            distFromRoot <- node.depth.edgelength(phy)
            phy$edge.length = (exp(param*distFromRoot[descendent])-exp(param*distFromRoot[parent]))/param
        }
    },
    "lambda"={
        # Pagel's lambda tree transformation
        if(param!=1) {
            if(is.ultrametric(phy)){
                root2tipDist <- max(branching.times(phy))
            }else{
                root2tipDist <- node.depth.edgelength(phy)[1:n] # for non-ultrametric trees. The 'up' limit should be exactly 1 to avoid singularity issues
            }
            phy$edge.length <- phy$edge.length * param
            phy$edge.length[extern] <- phy$edge.length[extern] + (root2tipDist * (1-param))
        }
        
    })
    
    # Add measurment error
    if(is.numeric(mserr)) phy$edge.length[extern] = phy$edge.length[extern] + mserr
    
    # Check for negative values
    if(any(phy$edge.length<0)) return(list(phy = phy, diagWeight = Inf, X=Inf, Y=Inf, Xw=Inf, Yw=Inf, det=Inf))
    
    # Compute the independent contrasts scores
    C <- pruning(phy, trans=FALSE)
    Xi <- crossprod(C$sqrtM, X)
    Yi <- crossprod(C$sqrtM, Y)
    
    # Return the determinant
    if(REML) deterM <- sum(log(C$varNode)) else deterM <- C$det
    
    # Return the score, variances, and scaled tree
    return(list(phy = phy, diagWeight = diagWeight, X=Xi, Y=Yi, Xw=X, Yw=Y, det=deterM))
}
