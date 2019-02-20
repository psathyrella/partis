################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################

# we can change to summary instead of print if we want to keep the same output between functions
print.fit_t.env<-function(x,...){
    cat("\n")
    if(is.function(x$model)){
    message("-- Summary results for user defined model --","\n")
    }else{
    message("-- Summary results for the ",x$model," model --","\n")
    }
    if(x$convergence==0){
        cat("successful convergence of the optimizer","\n")
    }else{
        cat("convergence of the optimizer has not been reached, try simpler model, different starting values or increase maxit","\n")
    }
    if(x$hess.value==0){
        cat("a reliable solution has been reached","\n")
    }else{
        cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")
    }
    cat("\n")
    cat("LogLikelihood:","\t",x$LH,"\n")
    cat("AIC:","\t",x$aic,"\n")
    cat("AICc:","\t",x$aicc,"\n")
    cat(x$free.parameters,"parameters","\n")
    cat("\n")
    if(!is.function(x$model)){
    cat("Beta:","\n")
    cat("______________________","\n")
    cat(x$param[2],"\n")
    cat("\n")
    cat("Brownian rate (sigma):","\n")
    cat("______________________","\n")
    cat(x$param[1],"\n")
    }else{
    cat("Parameters estimates:","\n")
    cat("______________________","\n")
    cat(x$param,"\n")
    }
    cat("\n")
    cat("Estimated root states:","\n")
    cat("______________________","\n")
    cat(x$root)
    cat("\n")
}