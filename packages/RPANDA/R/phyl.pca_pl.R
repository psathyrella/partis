################################################################################
##                                                                            ##
##                 RPANDA : (PL) ancestral states and PCA                     ##
##                                                                            ##
##   Julien Clavel - 01-02-2018                                               ##
##   S3 method for "ancestral"                                                ##
##   require: corpcor                                                         ##
################################################################################

# GLS estimation of ancestral states returned by the PL functions
ancestral.fit_pl.rpanda <- function(object){

  # extract objects
  if(!inherits(object,"fit_pl.rpanda")) stop("only works with \"fit_pl.rpanda\" class objects. See ?fit_t_pl")
  tree <- object$corrstruct$phy
  n <- object$n
  p <- object$p
  Y <- as.matrix(object$corrstruct$Yw)
  
  # covariance for the nodes
  V <- .vcvPhyloInternal(tree)
  indice <- (1:n)
  AY <- V[-indice,indice]
  vY <- V[indice,indice]

  # Ancestral state at the root
  one <- object$corrstruct$Xw
  a <- pseudoinverse(object$corrstruct$X)%*%object$corrstruct$Y
  
  # states at the nodes
  recons_t <- (AY%*%pseudoinverse(vY)%*%(Y-one%*%a))+(one[-1]%*%a)
  colnames(recons_t) = colnames(Y)
  rownames(recons_t) = paste("node_",n+1:Nnode(tree), sep="")

  # return the results
  res <- list(root=a, nodes=recons_t)
return(res)

}


## ------ Phylogenetic PCA on the penalized covariance matrix
# see also Revell 2009 - Evolution
phyl.pca_pl <- function(object, plot=TRUE, ...){
  
  if(!inherits(object,"fit_pl.rpanda")) stop("only works with \"fit_pl.rpanda\" class objects. See ?fit_t_pl")
  tree <- object$corrstruct$phy
  n <- object$n
  p <- object$p
  Y <- as.matrix(object$corrstruct$Yw)
  covR <- object$R$R
  
  # optional arguments
  args <- list(...)
  if(is.null(args[["axes"]])) axes <- c(1,2) else axes <- args$axes
  if(is.null(args[["col"]])) col <- "black" else col <- args$col
  if(is.null(args[["pch"]])) pch <- 19 else pch <- args$pch
  if(is.null(args[["cex"]])) cex <- 0.7 else cex <- args$cex
  if(is.null(args[["las"]])) las <- 1 else las <- args$las
  if(is.null(args[["main"]])) main <- "Regularized Phylogenetic PCA" else main <- args$main
  if(is.null(args[["mode"]])) mode <- "cov" else mode <- args$mode
  
  # if correlation matrix?
  if(mode=="corr") covR <- cov2cor(covR)
  
  # ancestral states
  anc <- ancestral(object)
  a <- anc$root
  nodes <- anc$nodes
  one <- object$corrstruct$Xw
  
  # PCA
  eig <- eigen(covR);
  U <- eig$vectors
  S <- (Y - one%*%a)%*%U # PCA scores for all tip species
  Srec <- (nodes - one[-1]%*%a)%*%U # PCA scores for the reconstructed states
  values <- eig$values
  
  # loadings (see Revell 2009)
  Cinv <- solve(vcv.phylo(tree))
  K <- (t(Y - one%*%a)%*%Cinv%*%S)/(n-1)
  L <- matrix(0,p,p);
  for(i in 1:p)
    for(j in 1:p)
      L[i,j]<-K[i,j]/sqrt(covR[i,i]*values[j])
  
  
  # plot
  if(plot){
    # contribution % variance
    tot<-sum(values)
    valX<-round(values[axes[1]]*100/tot,digits=2)
    valY<-round(values[axes[2]]*100/tot, digits=2)
    xlabel <- paste("PC",axes[1]," (",valX," %)", sep="")
    ylabel <- paste("PC",axes[2]," (",valY," %)", sep="")
    plot(S[,axes], main=main, xlab=xlabel, ylab=ylabel, pch=pch, col=col, las=las)
    abline(h=0,v=0)
    text(S[,axes],tree$tip.label, pos=2, cex=cex)
  }

  # results
  res <- list(values=values, scores=S, loadings=L, nodes_scores=Srec, mean=a, vectors=U)
  class(res) <- "pl.phyl.pca"
  invisible(res)
  #return(res) # don't need to return for plot only?
}
