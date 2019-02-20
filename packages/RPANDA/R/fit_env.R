fit_env <- function (phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par, df=NULL, f=1,
           meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE, expo.lamb=FALSE,
           expo.mu=FALSE, fix.mu=FALSE, dt=0, cond="crown")
{
  # first a spline is used to build the approximation model Env(t)
  if (is.null(df))
  {
    df <- smooth.spline(x=env_data[,1], env_data[,2])$df
  }
  spline_result <- sm.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  # In order to perform computation, the env_func is tabulated
  # control from lower_bound -10%, upper_bound + 10%
  lower_bound_control <- 0.10
  upper_bound_control <- 0.10
  lower_bound <- min(env_data[,1])
  upper_bound <- max(env_data[,1])
  # Tabulation of the function from lower_bound -10%, upper_bound + 10%
  time_tabulated <- seq(from=lower_bound*(1.0-lower_bound_control),
                        to=upper_bound*(1.0+upper_bound_control),
                        length.out=1+1e6)
  env_tabulated <- env_func(time_tabulated)
  # Tabulated function
  env_func_tab <- function(t)
  {
    b <- upper_bound * (1.0 + upper_bound_control)
    a <- lower_bound * (1.0 - lower_bound_control)
    # number of intervals
    n <- length(env_tabulated) - 1
    index <- 1 + as.integer( (t - a) * n / (b - a))
    return(env_tabulated[index])
  }
  f.lamb.env <- function(t,y){f.lamb(t, env_func_tab(t), y)}
  f.mu.env <- function(t,y){f.mu(t, env_func_tab(t), y)}
  res <- fit_bd(phylo, tot_time, f.lamb.env, f.mu.env, lamb_par, mu_par, f,
           meth, cst.lamb, cst.mu, expo.lamb, expo.mu, fix.mu, dt, cond)
  res$model <- "environmental birth death"
  res$f.lamb <- function(t){ abs(f.lamb(t, env_func_tab(t), res$lamb_par))}
  if (fix.mu==FALSE)
  {
    res$f.mu <- function(t){ abs(f.mu(t, env_func_tab(t), res$mu_par))}
  }
  class(res) <- "fit.env"
  return(res)
}
