##Univariate Discrete Distributions

|Functions|Emission Function Name|Parameters|
| ------- |-------------|----------|
|binomial_pdf (int k, int n, double p)|BINOMIAL|2|
|beta_binomial_pdf (int k, int n, double a, double b)|BETA_BINOMIAL|3|
|degenerate_pdf ( double value, double k)|DEGENERATE|1|
|discrete_uniform_pdf (int position, int a, int b)|DISCRETE_UNIFORM|2|
|hypergeometric_pdf (int k, int n, int N, int m)|HYPERGEOMETRIC|3|
|beta_negative_binomial_pdf (int k, int n, double a, double b)|BETA_NEGATIVE_BINOMIAL|3|
|maxwell_boltzman_pdf ( x, double a)|MAXWELL_BOLTZMAN|1|
|geometric_pdf (int k, double p)|GEOMETRIC|1|
|logarithmic_pdf (int k, double p)|LOGARITHMIC|1|
|negative_binomial_pdf (int k, int r, double p)|NEGATIVE_BINOMIAL|2|
|poisson_pdf (int k, double lambda)|POISSON|1|
|yule_simon_pdf (int k, double p)|YULE_SIMON|1|
|zipf_pdf (int k, int N, double s)|ZIPF|2|
|zipf_mandelbrot_pdf (int k, int N, double s, double q)|ZIPF-MANDELBROT|3|


##Univariate Continuous Distributions

|Functions|Emission Function Name|Parameters|
| ------- | ----------- | -------- |
|arcsine_pdf ( double x)|ARCSINE|0|
|beta_pdf (double x, double a, double b)|BETA|2|
|logit_normal_pdf (double x, double mu, double sigma)|LOGIT_NORMAL|2|
|continuous_uniform_pdf (double x, double a, double b)|CONTINUOUS_UNIFORM|2|
|kumaraswamy_pdf (double x, double a, double b)|KUMARASWAMY|2|
|raised_cosine_pdf (double x, double mu, double s)|RAISED_COSINE|2|
|triangular_pdf (double x, double a, double b, double c)|TRIANGULAR|3|
|truncated_normal_pdf (double x, double mu, double sd, double a, double b)|TRUNCATED_NORMAL|4|
|u_quadratic_pdf (double x, double a, double b)|U_QUADRATIC|2|
|wigner_semicircle_pdf (double x, double r)|WIGNER_SEMICIRCLE|1|
|beta_prime_pdf (double x, double a, double b)|BETA_PRIME|2|
|chi_pdf (double x, double k)|CHI|1|
|chi_squared_pdf (double x, double k)|CHI_SQUARED|1|
|inverse_chi_squared_pdf (double x, double v)|INVERSE_CHI_SQUARED|1|
|scaled_inverse_chi_squared_pdf (double x, double v, double sigma_sqrd)|SCALED_INVERSE_CHI_SQUARED|2|
|dagum_pdf (double x, double p, double a, double b)|DAGUM|3|
|exponential_pdf (double x, double lambda)|EXPONENTIAL|1|
|f_pdf (double x, double d1, double d2)|F_DIST|2|
|fishers_z_pdf (double x, double d1, double d2)|FISHERS_Z|2|
|folded_normal_pdf (double x, double mu, double sigma_sqrd)|FOLDED_NORMAL|3|
|frechet_pdf (double x, double alpha, double s, double m)|FRECHET|3|
|gamma_pdf (double x, double alpha, double beta)|GAMMA|2|
|inv_gamma_pdf (double x, double alpha, double beta)|INVERSE_GAMMA|2|
|half_normal_pdf (double x, double sigma)|HALF_NORMAL|1|
|inv_gaussian_pdf (double x, double mu, double lambda)|INVERSE_GAUSSIAN|2|
|levy_pdf (double x, double mu, double scale)|LEVY|2|
|log_cauchy_pdf (double x, double mu, double sigma)|LOG_CAUCHY|2|
|log_laplace_pdf (double x, double mu, double b)|LOG_LAPLACE|2|
|log_logistic_pdf (double x, double a, double b)|LOG_LOGISTIC|2|
|log_normal_pdf (double x, double mu, double sigma_sqrd)|LOG_NORMAL|2|
|pareto_pdf (double x, double alpha, double x_m)|PARETO|2|
|nakagami_pdf (double x, double mu, double w)|NAKAGAMI|2|
|rayleigh_pdf (double x, double sigma)|RAYLEIGH|1|
|gumbel_type_two_pdf (double x, double a, double b)|GUMBEL_TYPE_TWO|2|
|weibull_distribution (double x, double lambda, double k)|WEIBULL|2|
|cauchy_pdf (double x, double x_o, double gamma)|CAUCHY|2|
|gumbel_pdf (double x, double mu, double beta)|GUMBEL|2|
|generalized_normal_pdf (double x,  mu,  alpha,  beta)|GENERALIZED_NORMAL|3|
|hyperbolic_secant_pdf (double x)|HYPERBOLIC_SECANT|0|
|laplace_pdf (double x, double mu, double b)|LAPLACE|2|
|logistic_pdf (double x, double mu, double s)|LOGISTIC|2|
|standard_normal_pdf (double x)|STANDARD_NORMAL|0|
|normal_pdf (double x, double mu, double sigma)|NORMAL|2|
|students_t_pdf (double x, double v)|STUDENT_T|1|
|gumbel_type_one_pdf (double x, double a, double b)|GUMBEL_TYPE_ONE|2|
|generalized_extreme_value_pdf (double x, double mu, double sigma, double xi)|GENERALIZED_EXTREME_VALUE|3|
|generalized_pareto_pdf (double x, double mu, double sigma, double xi)|GENERALIZED_PARETO|3|


##Multivariate Continuous Distributions

|Functions|Emission Function Name|Parameters|
| ------- | ----------- | -------- |
|dirichlet_pdf (const std::vector<double>& x,const std::vector<double>&alpha)|DIRICHLET|Same as # of variables|
|multivariate_ewens_pdf(const std::vector<double>& x, const double theta)|EWENS|1|


#List of Probability Distribution Function (PDF.h and PDF.cpp)


###Univariate Discrete Distributions
```
double 	StochHMM::binomial_pdf (int k, int n, double p)
double 	StochHMM::beta_binomial_pdf (int k, int n, double a, double b)
double 	StochHMM::degenerate_pdf (double value, double k)
double 	StochHMM::discrete_uniform_pdf (int position, int a, int b)
double 	StochHMM::hypergeometric_pdf (int k, int n, int N, int m)
double 	StochHMM::poisson_binomial_pdf (int k, std::vector< double > &p)
double 	StochHMM::beta_negative_binomial_pdf (int k, int n, double a, double b)
double 	StochHMM::maxwell_boltzman_pdf (double x, double a)
double 	StochHMM::geometric_pdf (int k, double p)
double 	StochHMM::logarithmic_pdf (int k, double p)
double 	StochHMM::negative_binomial_pdf (int k, int r, double p)
double 	StochHMM::poisson_pdf (int k, double lambda)
double 	StochHMM::yule_simon_pdf (int k, double p)
double 	StochHMM::zipf_pdf (int k, int N, double s)
double 	StochHMM::zipf_mandelbrot_pdf (int k, int N, double s, double q)
```

###Univariate Continuous Distributions
```
double 	StochHMM::arcsine_pdf (double x)
double 	StochHMM::beta_pdf (double x, double a, double b)
double 	StochHMM::logit_normal_pdf (double x, double mu, double sigma)
double 	StochHMM::continuous_uniform_pdf (double x, double a, double b)
double 	StochHMM::kumaraswamy_pdf (double x, double a, double b)
double 	StochHMM::raised_cosine_pdf (double x, double mu, double s)
double 	StochHMM::triangular_pdf (double x, double a, double b, double c)
double 	StochHMM::truncated_normal_pdf (double x, double mu, double sd, double a, double b)
double 	StochHMM::u_quadratic_pdf (double x, double a, double b)
double 	StochHMM::wigner_semicircle_pdf (double x, double r)
double 	StochHMM::beta_prime_pdf (double x, double a, double b)
double 	StochHMM::chi_pdf (double x, double k)
double 	StochHMM::chi_squared_pdf (double x, double k)
double 	StochHMM::inverse_chi_squared_pdf (double x, double v)
double 	StochHMM::scaled_inverse_chi_squared_pdf (double x, double v, double sigma_sqrd)
double 	StochHMM::dagum_pdf (double x, double p, double a, double b)
double 	StochHMM::exponential_pdf (double x, double lambda)
double 	StochHMM::f_pdf (double x, double d1, double d2)
double 	StochHMM::fishers_z_pdf (double x, double d1, double d2)
double 	StochHMM::folded_normal_pdf (double x, double mu, double sigma_sqrd)
double 	StochHMM::frechet_pdf (double x, double alpha, double s, double m)
double 	StochHMM::gamma_pdf (double x, double alpha, double beta)
double 	StochHMM::inv_gamma_pdf (double x, double alpha, double beta)
double 	StochHMM::half_normal_pdf (double x, double sigma)
double 	StochHMM::inv_gaussian_pdf (double x, double mu, double lambda)
double 	StochHMM::levy_pdf (double x, double mu, double scale)
double 	StochHMM::log_cauchy_pdf (double x, double mu, double sigma)
double 	StochHMM::log_laplace_pdf (double x, double mu, double b)
double 	StochHMM::log_logistic_pdf (double x, double a, double b)
double 	StochHMM::log_normal_pdf (double x, double mu, double sigma_sqrd)
double 	StochHMM::pareto_pdf (double x, double alpha, double x_m)
double 	StochHMM::nakagami_pdf (double x, double mu, double w)
double 	StochHMM::rayleigh_pdf (double x, double sigma)
double 	StochHMM::gumbel_type_two_pdf (double x, double a, double b)
double 	StochHMM::weibull_distribution (double x, double lambda, double k)
double 	StochHMM::cauchy_pdf (double x, double x_o, double gamma)
double 	StochHMM::gumbel_pdf (double x, double mu, double beta)
double 	StochHMM::generalized_normal_pdf (double x, double mu, double alpha, double beta)
double 	StochHMM::hyperbolic_secant_pdf (double x)
double 	StochHMM::laplace_pdf (double x, double mu, double b)
double 	StochHMM::logistic_pdf (double x, double mu, double s)
double 	StochHMM::standard_normal_pdf (double x)
double 	StochHMM::normal_pdf (double x, double mu, double sigma)
double 	StochHMM::students_t_pdf (double x, double v)
double 	StochHMM::gumbel_type_one_pdf (double x, double a, double b)
double 	StochHMM::generalized_extreme_value_pdf (double x, double mu, double sigma, double xi)
double 	StochHMM::generalized_pareto_pdf (double x, double mu, double sigma, double xi)
```

###Multivariate Distributions
```
double 	StochHMM::dirichlet_pdf (const std::vector< double > &x, const std::vector< double > &alpha)
double 	StochHMM::multivariate_ewens_pdf (const std::vector< double > &x, double theta)
```

##Overloaded Functions for use with StochHMM::StateFuncs
###Univariate Discrete Distributions
```
double 	StochHMM::binomial_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::beta_binomial_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::degenerate_pdf (const double value, const std::vector< double > *param)
double 	StochHMM::discrete_uniform_pdf (const double position, const std::vector< double > *param)
double 	StochHMM::hypergeometric_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::beta_negative_binomial_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::maxwell_boltzman_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::geometric_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::logarithmic_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::negative_binomial_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::poisson_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::yule_simon_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::zipf_pdf (const double k, const std::vector< double > *param)
double 	StochHMM::zipf_mandelbrot_pdf (const double k, const std::vector< double > *param)
```
###Univariate Continuous Distributions
```
double 	StochHMM::arcsine_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::beta_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::logit_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::continuous_uniform_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::kumaraswamy_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::raised_cosine_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::triangular_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::truncated_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::u_quadratic_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::wigner_semicircle_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::beta_prime_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::chi_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::chi_squared_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::inverse_chi_squared_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::scaled_inverse_chi_squared_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::dagum_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::exponential_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::f_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::fishers_z_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::folded_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::frechet_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::gamma_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::inv_gamma_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::half_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::inv_gaussian_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::levy_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::log_cauchy_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::log_laplace_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::log_logistic_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::log_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::pareto_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::nakagami_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::rayleigh_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::gumbel_type_two_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::weibull_distribution (const double x, const std::vector< double > *param)
double 	StochHMM::cauchy_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::gumbel_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::generalized_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::hyperbolic_secant_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::laplace_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::logistic_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::standard_normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::normal_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::students_t_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::gumbel_type_one_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::generalized_extreme_value_pdf (const double x, const std::vector< double > *param)
double 	StochHMM::generalized_pareto_pdf (const double x, const std::vector< double > *param)
```

###Multivariate Continuous Distributions
```
double dirichlet_multi_pdf(const std::vector<double>* x, const std::vector<double>* param)
double ewens_multi_pdf(const std::vector<double>* x, const std::vector<double>* param)
```