plot_fit_bd <- function(fit.bd,tot_time)
{
  if (!inherits(fit.bd, "fit.bd"))
      stop("object \"fit.bd\" is not of class \"fit.bd\"")
  t <- seq(0,tot_time, length.out=100)
  dev.new()
  plot(-t, fit.bd$f.lamb(t), type='l', xlab="time", ylab="speciation rate", main="Fitted speciation rate")

  if ("f.mu" %in% attributes(fit.bd)$names)
  {
    # Attribute f.mu ==> not fixed extinction
    dev.new()
    plot(-t, fit.bd$f.mu(t), type='l', xlab="time", ylab="extinction rate", main="Fitted extinction rate")
    r <- function(t) {fit.bd$f.lamb(t) - fit.bd$f.mu(t)}
    dev.new()
    plot(-t, r(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
  }
  else
  {
    dev.new()
    plot(-t, fit.bd$f.lamb(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
  }
}
