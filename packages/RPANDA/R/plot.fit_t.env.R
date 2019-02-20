################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################

plot.fit_t.env<-function(x,steps=100,...){
    
    # Rates through time function
    if(is.function(x$model)){
        fun_temp<-function(x,temp,model,param){
            rate_fun<-function(x){ model(x, temp, param) }
            rate<-rate_fun(x)
            return(rate)
        }
    }else if(x$model=="EnvExp"){
        fun_temp<-function(x,temp,model,param){
            sig<-param[1]
            beta<-param[2]
            rate<-(sig*exp(beta*temp(x)))
            return(rate)
        }
    }else if(x$model=="EnvLin"){
        fun_temp<-function(x,temp,model,param){
            sig<-param[1]
            beta<-param[2]
            rate<-sig+(beta-sig)*temp(x)
            return(rate)
        }
    }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=steps)
    
    # Rates through time
    rates <- fun_temp( x=t, temp=x$env_func, model=x$model, param=x$param)
    
    plot(-t, rates, type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), ...)
    results<-list(time_steps=t, rates=rates)
    invisible(results)
}

# Allows drawing lines and superposing various results

lines.fit_t.env<-function(x,steps=100,...){
    
    # Rates through time function
    if(is.function(x$model)){
        fun_temp<-function(x,temp,model,param){
            rate_fun<-function(x){ model(x, temp, param) }
            rate<-rate_fun(x)
            return(rate)
        }
    }else if(x$model=="EnvExp"){
        fun_temp<-function(x,temp,model,param){
            sig<-param[1]
            beta<-param[2]
            rate<-(sig*exp(beta*temp(x)))
            return(rate)
        }
    }else if(x$model=="EnvLin"){
        fun_temp<-function(x,temp,model,param){
            sig<-param[1]
            beta<-param[2]
            rate<-sig+(beta-sig)*temp(x)
            return(rate)
        }
    }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=steps)
    
    # Rates through time
    rates<-fun_temp( x=t, temp=x$env_func, model=x$model, param=x$param)
    
    lines(-t, rates, type='l', ...)
    results<-list(time_steps=t, rates=rates)
    invisible(results)
}