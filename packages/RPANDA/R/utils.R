################################################################################
##                                                                            ##
##                                RPANDA : Utils                              ##
##                                                                            ##
##   Julien Clavel - 01-02-2018                                               ##
##   S3 methods, simulations, miscellaneous                                   ##
##                                                                            ##
################################################################################

# S3 generic method "ancestral" for reconstructing or retrieving ancestral states (see phyl.pca_pl.R) # should I use "predict" instead?
ancestral <- function(object) UseMethod("ancestral")

# S3 for the fit_t_env class
ancestral.fit_t.env <- function(object){
    
    # extract objects
    if(!inherits(object,"fit_t.env")) stop("only works with \"fit_t.env\" class objects. See ?fit_t_env")
    
    # Ancestral state at the root
    a <- object$root
    names(a) = "root"

    res <- a
 return(res)
 warning("only the root state is currently estimated for models of the class \"fit_t.env\"") # To remove later

}

# S3 for the fit_t_comp class?
# TODO
ancestral.fit_t.comp <- function(object){
    
    # extract objects
    if(!inherits(object,"fit_t.comp")) stop("only works with \"fit_t.comp\" class objects. See ?fit_t_comp")
    
    anc <- object$z0
    names(anc) ="root"
    return(anc)
    warning("only the root state is currently estimated for models of the class \"fit_t.comp\"")
}

# Build a matrix with tip and internal covariances
.vcvPhyloInternal <- function(tree){
    nbtip <- Ntip(tree)
    dis <- dist.nodes(tree)
    MRCA <- mrca(tree, full = TRUE)
    M <- dis[as.character(nbtip + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    return(M)
}

# Build the matrix square root and inverse matrix square root of the phylogenetic covariance matrix
.transformsqrt <- function(tree){
    vcv_tr <- vcv.phylo(tree)
    eig <- eigen(vcv_tr)
    sqrtmValues <- sqrt(eig$values)
    sqrtM1 <- tcrossprod(eig$vectors%*%diag(1/sqrtmValues), eig$vectors)
    sqrtM <- tcrossprod(eig$vectors%*%diag(sqrtmValues), eig$vectors)
    matrices <- list(sqrtM1=sqrtM1, sqrtM=sqrtM)
    return(matrices)
}

# ---- Function to simulate random covariance matrices with a specified eigenstructure
# From Uyeda et al. 2015 - Systematic Biology 64(4):677-689.
Posdef <- function (p, ev = rexp(p, 1/100)) {
  Z <- matrix(ncol=p, rnorm(p^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
