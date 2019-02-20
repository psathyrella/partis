plot_spectR <- function(spectR)
{
  if (!inherits(spectR, "spectR"))
      stop("object \"spectR\" is not of class \"spectR\"")
  
  m=subset(spectR$eigenvalues,spectR$eigenvalues>=1)	
		#get spectral density
		dens(m)->d
			integr(d$x,d$y)->dint
				(d$y/dint)->dsc

			
		par(mfrow=c(1,2))
			plot(d$x,dsc,type="l",ann=F)
				mtext(expression(f*(x)/integral(f(y)*dy)),2,2)
					mtext("ln eigenvalue",1,2)
			plot(sort(log(m),decreasing=T),ann=F)
				mtext("rank",1,2)
					mtext("ln eigenvalue",2,3)				
}

#gaussian kernel convolution		
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    na.omit(x)->x
    if(length(x) == 0)
        stop('too infinite.')
  }
  	kernelG<-function(x, mean=0, sd=1) 
		dnorm(x, mean = mean, sd = sd)
	x <- log(x)	
	sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
	X <- seq(from, to, len = n)
	M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

#integration
integr <- function(x, f)
{
       if (!is.numeric(x))
       {
              stop('"x" is not numeric.')
       }
       if (!is.numeric(f))
       {
              stop('"f" is not numeric.')
       }
       if (length(x) != length(f))
       {
              stop('integration variable and integrand are wrong for each other.')
       }

       length(x)->n
       integral=0.5*sum((x[2:n]-x[1:(n-1)])*(f[2:n]+f[1:(n-1)]))
       return(integral)
}



