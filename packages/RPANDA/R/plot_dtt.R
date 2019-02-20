plot_dtt <- function(fit.bd,tot_time,N0)
{
  if (!inherits(fit.bd, "fit.bd"))
      stop("object \"fit.bd\" is not of class \"fit.bd\"")
  t <- seq(0,tot_time, length.out=100)

  if ("f.mu" %in% attributes(fit.bd)$names)
  {
    # Attribute f.mu ==> not fixed extinction
    r <- function(t) {-fit.bd$f.lamb(t) + fit.bd$f.mu(t)}
    R <- function(s){.Integrate(Vectorize(r),0,s)}
    N <- N0 * exp(Vectorize(R)(t))
    dev.new()
    plot(-t, N, type='l', xlab="time", ylab="Number of species", main="Diversity Through Time")
  }
  else
  {
    r <- function(t) {-fit.bd$f.lamb(t)}
    R <- function(s){.Integrate(Vectorize(r),0,s)}
    N <- N0 * exp(Vectorize(R)(t))
    dev.new()
    plot(-t, N, type='l', xlab="time", ylab="Number of species", main="Diversity Through Time")
  }
}
