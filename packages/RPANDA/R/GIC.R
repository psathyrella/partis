################################################################################
##                                                                            ##
##  RPANDA : Generalized Information Criterion (GIC) for high-dimensional     ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################



gic_criterion <- function(Y, tree, model="BM", method=c("RidgeAlt","RidgeArch","LASSO","ML","RidgeAltapprox","LASSOapprox"), targM=c("null","Variance","unitVariance"), param=NULL, tuning=0, REML=TRUE, ...){
  
  # ellipsis for additional arguments
  par <- list(...)
  if(is.null(par[["corrstruct"]])){ corrstruct <- NULL }else{ corrstruct <- par$corrstruct}
  if(is.null(par[["SE"]])){ SE <- NULL}else{ SE <- par$SE}
  if(is.null(par[["scale.height"]])){ scale.height <- FALSE }else{ scale.height <- par$scale.height}
  
  # Scale the tree to unit length
  if(scale.height==TRUE) tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
  
  # Select the method
  method <- match.arg(method)[1]
  if(method=="RidgeAltapprox"){
      method <- "RidgeAlt"
  }else if(method=="LASSOapprox"){
      method <- "LASSO"
  }
  
  # Select the target
  targM <- match.arg(targM)[1]
  
  # Parameters
  n <- nrow(Y)
  p <- ncol(Y)
  
  # Y must be a matrix (force coercion)
  Y <- as.matrix(Y)
  
  # check for parameters
  if(is.null(param) & model!="BM") stop("please provide a parameter value for the evolutionary model!!")
  if(method=="ML" & p>=n) warning("The covariance matrix is singular, the log-likelihood (and the GIC) is unreliable!!")
  
  # Choose the model
  if(model=="BM") mod.par = 0 else mod.par = 1
  if(!is.null(SE)) mod.par = mod.par + 1
  nC <- n
  if(REML==TRUE) n <- n-1
  
  # Transform the tree if it's not a model fit
  if(is.null(corrstruct)){
      corrstruct <- .transformTree(tree, param, model=model, mserr=SE, Y=Y, X=matrix(1,nrow(Y)), REML=REML)
  }
  
  # square root transform
  D <- .transformsqrt(corrstruct$phy)$sqrtM1
  X <- crossprod(D, corrstruct$Xw)
  Y <- crossprod(D, corrstruct$Yw)
  beta <- pseudoinverse(X)%*%Y
  
  # Covariance matrix
  Yk <- Y - X%*%beta
  S <- crossprod(Yk)/n
  
  # Determinant for the phylogenetic tree
  Ccov <- corrstruct$det
  
  # Switch between targets
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
           P <- .makePenalty(S,tuning,Target,targM)$S
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
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
  
  # GIC score
  
  if(method=="RidgeArch"){
      
    # First and second derivative of the functional (we can use patterned matrix to target some matrix elements)
    # We use the Kronecker-vec identity to speed up the computations
    T1 <- sapply(1:nC, function(i){
        Sk <- tcrossprod(Yk[i,]) ;
        VSV <- 0.5*(P - (1-tuning)*Sk - tuning*Target);
        VSV2 <- 0.5*(P - Sk);
        sum(VSV * 2*(Pi%*%VSV2%*%Pi))
    })
    
    df = sum(T1)/nC
    sigma_df <- df
    
  }else if(method=="LASSO" | method=="ML"){
    # LASSO or ML
    Tf2 <- function(S, P) {
      I <- ifelse(P==0,0,1) ;
      t(.vec(S*I))%*%.vec(P%*%(S*I)%*%P)
    }
    
    sigma_df <- (1/(2*nC))*sum(sapply(1:nC, function(i){
      Sk <- Yk[i,]%*%t(Yk[i,]) ;
      Tf2(Sk, Pi)})) - (1/2)*Tf2(S,Pi)
    
  }else if(method=="RidgeAlt"){
    # Alternative Ridge
    H <- (1/(0.5*(kronecker(d,d)+tuning)))
    
    # 2) First derivative of the functional
    T1 <- sapply(1:nC, function(i){
       Sk <- tcrossprod(Yk[i,]) ;
       VSV <- .vec(crossprod(V, (0.5*(P - (Sk - tuning*Target) - tuning*Pi))%*%V));
       VSV2 <- .vec(crossprod(V, (0.5*(P - Sk))%*%V));
      sum(VSV * (H*VSV2))
    })
    
    df = sum(T1)/nC
    sigma_df <- df
  }
  
  # Number of parameters for the root state:
  # The Information matrix from the Hessian and gradients scores
  XtX <- solve(crossprod(X))
  T2 <- sapply(1:nC, function(i){
      gradient <- X[i,]%*%t(Pi%*%t(Y[i,] - X[i,]%*%beta))
      sum(gradient * (XtX%*%gradient%*%P))
  })
  beta_df <- sum(T2)
  
  # LogLikelihood (minus)
  DP <- as.numeric(determinant(P)$modulus)
  llik <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*DP + n*sum(S*Pi))
  GIC <- 2*llik + 2*(sigma_df+beta_df+mod.par)
  
  # return the results
  results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=nC, bias=sigma_df+beta_df+mod.par)
  class(results) <- "gic.rpanda"
  return(results)
}

# Extractor for fit_pl.rpanda 'class'
GIC.fit_pl.rpanda <- function(object, ...){
  if(!inherits(object,"fit_pl.rpanda")) stop("only works with \"fit_pl.rpanda\" class objects")
  gic_criterion(Y=object$variables$Y, tree=object$variables$tree, model=object$model, method=object$method, targM=object$targM, param=object$model.par, tuning=object$gamma, REML=object$REML, corrstruct=object$corrstruct, SE=object$SE)
}
