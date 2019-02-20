################################################################################
##                                                                            ##
##  RPANDA : Print functions for the GIC and Penalized Likelihood models      ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################

# we can change to summary instead of print if we want to keep the same output between functions
print.fit_pl.rpanda<-function(x,...){
  cat("\n")
    message("-- Summary results for the ",x$model," model --","\n")
    message("Penalization: ",x$method,"\n")
  cat("LOOCV (negative):","\t",x$loocv,"\n")
  cat("\n")
    cat("Model parameter:","\n")
    cat("______________________","\n")
    cat(x$model.par,"\n")
    cat("\n")
    cat("Regularization parameter (gamma):","\n")
    cat("______________________","\n")
    cat(x$gamma,"\n")
    cat("\n")
    cat("Evolutionary Covariance of size:",x$p,"by",x$p,"\n")
    cat("for",x$n,"species","\n")
    cat("______________________","\n")
  cat("\n")
}

# GIC printing options
print.gic.rpanda<-function(x,...){
  cat("\n")
  message("-- Generalized Information Criterion --","\n")
  cat("GIC:",x$GIC,"| Log-likelihood",x$LogLikelihood,"\n")
  cat("\n")
}