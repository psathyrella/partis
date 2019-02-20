plot_fit_env <- function(fit.env, env_data, tot_time)
{
  if (!inherits(fit.env, "fit.env"))
      stop("object is not of class \"fit.env\"")
  t <- seq(0,tot_time, length.out=100)
  dev.new()
  plot(-t, fit.env$f.lamb(t), type='l', xlab="time", ylab="speciation rate", main="Fitted speciation rate")
  # Plot f.lamb(env_data)
  df <- smooth.spline(env_data[,1], env_data[,2])$df
  spline_result <- sm.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  dev.new()
  plot(env_func(t), fit.env$f.lamb(t), xlab="Environmental data", ylab="speciation rate", main="Fitted speciation rate")


  if ("f.mu" %in% attributes(fit.env)$names)
  {
    # Attribute f.mu ==> not fixed extinction
    dev.new()
    plot(-t, fit.env$f.mu(t), type='l', xlab="time", ylab="extinction rate", main="Fitted extinction rate")
    dev.new()
    plot(env_func(t), fit.env$f.mu(t), xlab="Environmental data", ylab="extinction rate", main="Fitted extinction rate")
    r <- function(t) {fit.env$f.lamb(t) - fit.env$f.mu(t)}
    dev.new()
    plot(-t, r(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
  dev.new()
    plot(env_func(t), r(t), xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
  else
  {
    dev.new()
    plot(-t, fit.env$f.lamb(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
    dev.new()
    plot(env_func(t), fit.env$f.lamb(t), xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
}
