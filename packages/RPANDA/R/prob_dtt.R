prob_dtt<- function(fit.bd, tot_time, time, N0, l=N0, f = l/N0, m = seq(N0), method="simple", lin = FALSE, prec = 1000,type = "stem",logged = TRUE){
  
  time.forward <- abs(time[order(abs(time))])
  
  tryIntegrate <- function(...){
    res <- try(integrate(...));
    if (class(res) == 'try-error'){
      return(-Inf)
    }else{
      return(res$value)
    }
  }
  
  qu.gen <- function(a,b,f.lamb,f.mu,method,prec){
    r<-function(u){f.lamb(u)-f.mu(u)} 
    r.int<-function(x,y){tryIntegrate(Vectorize(r),x,y,stop.on.error=FALSE)} 
    r.int.q<-function(y){exp(-r.int(a,y))*f.mu(y)} 
    r.int.int<-function(x,y){tryIntegrate(Vectorize(r.int.q),x,y,stop.on.error=FALSE)}
    int <- r.int.int(a,b) 
    if (method!="simple") int <- mpfr(int,prec)/(1+mpfr(int,prec))     
    return(int)
  }
  
  eta.gen <- function(a,b,f.lamb,f.mu,method,prec){
    r<-function(x){f.lamb(x)-f.mu(x)} 
    r.int<-function(x,y){tryIntegrate(Vectorize(r),x,y,stop.on.error=FALSE)} 
    r.int.n<-function(y){exp(r.int(y,b))*f.lamb(y)} 
    r.int.int<-function(x,y){tryIntegrate(Vectorize(r.int.n),x,y,stop.on.error=FALSE)}
    int <- r.int.int(a,b)
    if (method!="simple") int <- mpfr(int,prec)/(1+mpfr(int,prec))     
    return(int)
  }
  
  stem.prob <- function(y, qst, nst, logged, method){
    if (method=="simple"){
      prob <- log(1/(1+qst)) + log(1/(1+nst)) + (y-1)*log(nst/(1+nst))
    }else{        
      prob <- (1-qst)*(1-nst)*(nst)^(y-1)
    }
    if(!logged & method=="simple") prob <- exp(prob)
    return(prob)
  }
  
  crown.prob <- function(y, qst, nst, logged, method){
    if (method=="simple"){
      sm <- y-1 + (2 * (qst/(1+qst) * nst/(1+nst))) / ((1/(1+qst)) * (1/(1+nst)))
      prob <- 2*log(1/(1+qst)) + 2*log(1/(1+nst)) + (y-2)*log(nst/(1+nst)) + log(sm)
    }else{
      sm <- (y-1) + (2 * qst * nst) / ((1-qst) * (1-nst))
      prob <- (1-qst)^2*(1-nst)^2*(nst)^(y-2)*sm
    }        
    if(!logged & method=="simple") prob <- exp(prob)
    return(prob)
  }
  
  Init <- function(m,type){
    p <- mat.or.vec(nr = length(m), nc= 1)
    if (type=="stem" & m[1]==1) p[1] <- 1 
    if (type=="crown" & m[1]==2) p[1] <- 1 
    return(p)
  }
  
  Final <- function(m,N0){
    if (length(m)==1){
      p <- ifelse(m==N0,1,0) 
    }else{
      p <- mat.or.vec(nr = length(m), nc= 1)
      if (N0<=max(m)) p[N0] <- 1 
    }
    return(p)
  }
  
  Prob.t.allsp <- function(n, m, s=0, t, tot_time, Pnx, f.lamb, f.mu, type="stem", logged=TRUE, method){
    
    Pnm.prob <- function(n, x, qst, nst, logged){ 
      Pnm.prob.1 <- function(n,x,qst,nst,logged){
        if(x > n) prob <- -Inf*nst else prob <- x * log(1/(1+nst)) + (n-x) *log(nst/(1+nst)) + lchoose(n-1, x-1)
        if (!logged) prob <- exp(prob)
        return(prob)
      }
      
      Pnm.prob.2 <- function(n,x,qst,nst,logged){
        left  <- x * log(1/(1+qst)) + x * log(1/(1+nst)) + (n-x) * log(nst/(1+nst)) + lgamma(x+1) + lgamma(n)
        base  <- log(qst/(1+qst)) + log(nst/(1+nst)) - log(1/(1+nst)) - log(1/(1+qst))
        ifelse(x <= n, k <- seq(0, x - 1), k <- seq(x - n, x - 1)) 
        
        Pnm.prob.2.1 <- function(base,k,n,x) return(k * base - lgamma(k+1)-lgamma(x-k+1)-lgamma(x-k)-lgamma(n+k-x+1))
        
        Pnm.prob.2.2 <- function(pk, left, logged){
          if(!is.matrix(pk)){ 
            tot <- exp(pk) 
            tot <- ifelse(tot==Inf, tot <- as.numeric(exp(mpfr(pk,1000))),tot)
            prob <- left + log(tot)
          }else{
            tot <- colSums(exp(t(t(pk)-colMeans(pk))))
            tot <- ifelse(tot==Inf, as.numeric(colSums(exp(t(t(pk)-mpfr(colMeans(pk),1000))))),tot) 
            prob <- left + log(tot)+colMeans(pk)
          }
          if (!logged) prob <- exp(prob)
          return(prob)
        }
        
        pk <- mapply(Pnm.prob.2.1, base, MoreArgs = list(k,n,x))
        prob <- Pnm.prob.2.2(pk=pk, left=left, logged=logged)
      }
      
      prob <- ifelse((qst == 0), Pnm.prob.1(n,x,qst,nst,logged), Pnm.prob.2(n,x,qst,nst,logged))      
      return(prob) # Return probability
    }
    
    p <- mat.or.vec(nr=length(m), nc=length(t))
    init.final <- t==0|t==tot_time
    if (init.final[1]==T) p[,1] <- Init(m, type=type)
    if (init.final[length(t)]==T) p[,length(t)] <- Final(m = m, N0 = N0)      
    t <- t[!init.final]
    
    qst.tp <- mapply(qu.gen,a=t,MoreArgs= list(b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec))
    nst.tp <- mapply(eta.gen, a=t, MoreArgs=list(b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec)) 
    qst.st  <- mapply(qu.gen, b=t, MoreArgs= list(a=s,f.lamb=f.lamb, f.mu=f.mu,method=method,prec=prec))
    nst.st 	<- mapply(eta.gen, b=t, MoreArgs = list(a=s, f.lamb=f.lamb, f.mu=f.mu,method=method,prec=prec))
    Pnm <- mapply(Pnm.prob, x=m, MoreArgs = list(n=n,qst=qst.tp, nst=nst.tp, logged=logged))
    if (type == "crown"){
      Pmx <- mapply(crown.prob, y=m, MoreArgs = list(qst=qst.st, nst=nst.st, logged=logged,method=method))
    } else {
      if (type == "stem"){
        Pmx <- mapply(stem.prob, y=m, MoreArgs= list(qst=qst.st, nst=nst.st, logged=logged, method=method))
      } else { stop("Needs to be either stem or crown")
      }
    }
    if(logged){
      mat <- exp(Pnm + Pmx - Pnx)
    } else {
      mat	<- (Pnm * Pmx) / Pnx
    }
    if (sum(init.final)==2){
      p[,2:(length(init.final)-1)] <- t(mat)
    }else if (sum(init.final)==1){
      if (init.final[1]==T) p[,2:length(init.final)] <- t(mat)
      if (init.final[length(t)]==T) p[,1:(length(init.final)-1)] <- t(mat)
    }else{
      p <- t(mat)
    }
    return(p)
  }
  
  F.calc <- function(z, l, m, qst, nst, F1=NULL, pb=NULL){
    prob.part.ext <- function(z, l, m, qst, nst){
      top <- (qst + z * (1-nst - qst))^m * (nst)^l
      bot <- (1 - z * nst)^(l + m)
      out <- top / bot
      return(out)
    }
    prob.part.sum <- function(z, j, qst, nst){
      out <- (((1- nst - qst) * (1 - z * nst))/((qst + z * (1-nst-qst)) * nst))^j
      return(out)
    }
    fact.part <- function(m, l, j){
      k 	<- l - j
      out <-chooseMpfr(m,j)*chooseMpfr(m+k-1,k)  
      return(out)
    }
    prob.ext <- prob.part.ext(z=z, l=l, m=m, qst=qst, nst=nst)
    mxsum <- min(m,l) 
    tot <-0
    y <- 0:mxsum
    if (is.null(pb) | is.null(F1)){
      tot <- sum(fact.part(m=m,l=l,j=y)*prob.part.sum(z=z, j=y, qst=qst, nst=nst))
      return(tot*prob.ext)
    }else{
      fact.part1 <- chooseMpfr(m,y)
      fact.part2 <- chooseMpfr(m+l-y-1,l-y)
      tot <- sum(pmax(fact.part1,fact.part2)/F1*pb*prob.part.sum(z=z, j=y, qst=qst, nst=nst)*prob.ext*pmin(fact.part1,fact.part2)) 
      return(tot)
    }
  }
  
  Prob.t.f <- function(f, l, m, s=0, t, tot_time, z, F1, f.lamb, f.mu, type, prec){      
    if(t==0){
      p <- Init(m,type)
    }else if (t==tot_time){
      p <- Final(m,N0)
    }else{
      qst.tp <- qu.gen(a=t,b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec)
        nst.tp <- eta.gen(a=t,b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec) 
        qst.st <- qu.gen(a=s, b=t, f.lamb=f.lamb, f.mu=f.mu,method=method,prec=prec)
        nst.st  <- eta.gen(a=s, b=t, f.lamb=f.lamb, f.mu=f.mu,method=method,prec=prec)
        if (type == "crown"){
          pb <- crown.prob (y=m[1], qst=qst.st, nst=nst.st, logged=logged, method=method)
          if (length(m)>1) pb <- c(pb,sapply(m[-1],function(y) crown.prob (y=y, qst=qst.st, nst=nst.st, logged=logged, method=method)))
        } else{
          pb <- stem.prob (y=m[1], qst=qst.st, nst=nst.st, logged=logged,method=method)
          if (length(m)>1) pb <- c(pb,sapply(m[-1],function(y) stem.prob (y=y, qst=qst.st, nst=nst.st, logged=logged, method=method)))
        }
        p <- as.numeric(F.calc(z=z, l=l, m=m[1], qst=qst.tp, nst=nst.tp, pb=pb[1], F1=F1))
        if (length(m)>1) p <- c(p,sapply(m[-1],function(y) as.numeric(F.calc(z=z, l=l, m=y, qst=qst.tp, nst=nst.tp, pb=pb[which(m==y)], F1=F1))))
      }
      return(p)
    }  
    
    if (!(method=="simple"|method=="hard")) stop ("method needs to be either \'simple\' or \'hard\'")
    
    if (!inherits(fit.bd, "fit.bd")) 
      stop("object \"fit.bd\" is not of class \"fit.bd\"")
    if ("f.mu" %in% attributes(fit.bd)$names==FALSE) {
      f.mu <- function(t){0}
    }else{
      f.mu <- function(t) abs(fit.bd$f.mu(tot_time-t))}
      f.lamb <- function(t) abs(fit.bd$f.lamb(tot_time-t))
    
    if(type == "crown"){
      m <- m[m > 1] 
      if(length(m) == 0) stop('Probability undefined for less than 2 species') 
    }
    s=0
    qst.sp <- qu.gen(a=s,b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec)
    nst.sp <- eta.gen(a=s,b=tot_time,f.lamb=f.lamb,f.mu=f.mu,method=method,prec=prec) 
    z <- 1 - f
    if (method=="simple"){
      if (type == "crown"){
        Pnx <- crown.prob(y=N0, qst=qst.sp, nst=nst.sp, logged=logged,method=method)
      } else {
        if (type == "stem"){
          Pnx <- stem.prob(y=N0, qst=qst.sp, nst=nst.sp, logged=logged,method=method)
        } else { stop("Needs to be either stem or crown")
        }
      }
      Prob <-  Prob.t.allsp(n=N0, m=m, s=0, t=time.forward, tot_time=tot_time, Pnx=Pnx, f.lamb=f.lamb, f.mu=f.mu, type=type, logged=logged, method=method)
    }else{
      if (type == "crown") {
        F1 <- F.calc(z=z, l=l, m=2, qst=qst.sp, nst=nst.sp, F1=NULL, pb=NULL)
      } else {
        if (type != "stem") stop("'type' needs to be specified as either stem or crown")
        F1 <- F.calc(z=z, l=l, m=1, qst=qst.sp, nst=nst.sp, F1=NULL, pb=NULL)
      }
      Prob <- matrix(nrow = length(m), ncol=length(time.forward))
      for (i in 1:length(time.forward)){
        if (length(m)==1){
          Prob[i]<- Prob.t.f(f=f, l=l, m=m, s=0, t=time.forward[i], tot_time=tot_time, z=z, F1=F1, f.lamb=f.lamb, f.mu=f.mu, type=type, prec = prec)
        }else{
          Prob[,i]<- Prob.t.f(f=f, l=l, m=m, s=0, t=time.forward[i], tot_time=tot_time, z=z, F1=F1, f.lamb=f.lamb, f.mu=f.mu, type=type, prec = prec)
        }
      }
    }
    if(length(m)>1){
      rownames(Prob) <- m
      colnames(Prob) <- tot_time-time.forward
    }
    return(Prob)
  }
  
