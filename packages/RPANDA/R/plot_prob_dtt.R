plot_prob_dtt <- function (mat, grain=0.1, plot.prob=TRUE, plot.mean = TRUE, int = TRUE, plot.bound=FALSE, conf = 0.95, add=FALSE, col.mean="red", col.bound="blue", lty="solid", lwd=1){
  int.conf <- function (vect, conf = 0.95, div){
    tab <- cbind(vect,div)
    tab.sort <- tab[order(tab[,1], decreasing = TRUE),]
    test = tab.sort[1,1]
    i = 1
    while (test < conf & i+1 <= nrow(tab.sort)){
      test = test+tab.sort[i+1,1]
      i = i+1
    }
    out <- ifelse(tab[,2]<min(tab.sort[c(1:i),2])|tab[,2]>max(tab.sort[c(1:i),2]),0,tab[,1])
    return(out)
  }
  
  nsp <- nrow(mat)
  tim <- ncol(mat)
  pm   <- (matrix(nrow = nsp * tim, ncol = 4))
  pm[,1] <- as.numeric(rownames(mat)) 
  pm[,2] <- sapply(as.numeric(colnames(mat)), rep, times=nsp)
  if (int){
    mat.int.conf <- matrix(nrow=nsp, ncol = tim)
    mat.int.conf <- apply(mat,2,function(y) int.conf(y,conf=0.95, div = as.numeric(rownames(mat))))
    pm[,3] <- mat.int.conf
    pm[,4] <- c('white',paste0('gray',99:1),'black')[as.numeric(cut(mat.int.conf,breaks=c(seq(0,grain,length.out=99),0.5,1)))]
    if(plot.bound){
      bound.min <- pm[apply(X = mat.int.conf, MARGIN = 2,function(y) min(which(y!=0))),1]
      bound.max <- pm[apply(X = mat.int.conf, MARGIN = 2,function(y) max(which(y!=0))),1]
    }
  }else{
    pm[,3] <- mat
    pm[,4] <- c('white',paste0('gray',99:1),'black')[as.numeric(cut(mat,breaks=c(seq(0,grain,length.out=99),0.5,1)))]    
  }
  if (add==TRUE){
    if (plot.prob==TRUE) points(pm[,2],pm[,1],col=pm[,4],pch=20)
  }else{
    if (plot.prob==TRUE) plot(pm[,2],pm[,1],col=pm[,4],pch=20, xlim=c(max(as.numeric(pm[,2])),min(as.numeric(pm[,2]))),ann=FALSE)
  }
  if(plot.bound){
    if(plot.prob==FALSE){
      if(add==FALSE){
        plot(as.numeric(colnames(mat)), bound.max,type = "l", col = col.bound,ann=FALSE, lwd=lwd,xlim=c(max(as.numeric(pm[,2])),min(as.numeric(pm[,2]))))
        lines(as.numeric(colnames(mat)), bound.min,col = col.bound, lwd=lwd)
      }
      else{
        lines(as.numeric(colnames(mat)), bound.max,col = col.bound, lwd=lwd)
        lines(as.numeric(colnames(mat)), bound.min,col = col.bound, lwd=lwd)
      }
    }else{
      lines(as.numeric(colnames(mat)), bound.max,col = col.bound, lwd=lwd)
      lines(as.numeric(colnames(mat)), bound.min,col = col.bound, lwd=lwd)
    }
  }   
  if (plot.mean){
    mean.val <- t(as.numeric(rownames(mat))) %*% mat
    if (plot.prob==FALSE & plot.bound==FALSE){
      if (add==FALSE){
        plot(as.numeric(colnames(mat)), mean.val,type = "l", col = col.mean,ann=FALSE, lwd=lwd,xlim=c(max(as.numeric(pm[,2])),min(as.numeric(pm[,2]))))
      }else{
        lines(as.numeric(colnames(mat)), mean.val,col = col.mean, lty=lty, lwd=lwd)
      }
    }else{
      lines(as.numeric(colnames(mat)), mean.val,col = col.mean, lty=lty, lwd=lwd)
    }
  }
}
