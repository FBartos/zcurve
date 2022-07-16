/*ZcurveEM
* author: Frantisek Bartos
* email:  f.bartos96@gmail.com
*
* functions to fit z-curve using EM algorithm
* the most computation demanding part is written in Rcpp
*
* parameters:
* z:
*   z-values imput
* type:
*   1 = components with fixed means
*   2 = components with estimated means
* K:
*   number of components in case of components with estimaten means
* mu:
*   location of component's mean in case of their fixation
* theta_alpha:
*   alpha parameter of dirichlet distibution used to generate initial weights in the z-curve initialization process
* mu_alpha:
*   alpha parameter of dirichlet distibution used to generate initial means in the z-curve initialization process (only if type = 2)
* mu_max:
*   maximum value for mean component in the z-curve initialization process (only if type = 2)
* sig_level:
*   level of significance used, the lower censoring
* b:
*   the upper limit to which the z-curve is fitted (upper censoring)
* criterion:
*   criterion to terminate the z-curve fitting process
* max_iter:
*   maximum number of iteration in the z-curve fitting process
* criterion_start
*   criterion to terminate the z-curve initialization process
* max_iter_start:
*   maximum number of iteration in the z-curve initialization process
* fit_reps:
*   number of fits in the z-curve initialization process
* bootstrap:
*   number of bootstrap samples for confidence interval computation (FALSE for no bootstrap)
*/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.zdist_lpdf)]]
NumericVector zdist_lpdf(NumericVector x, double mu, double sigma, double a, double b) {
  NumericVector l1    = Rcpp::dnorm( x, mu, sigma, false);
  NumericVector l2    = Rcpp::dnorm(-x, mu, sigma, false);

  double l1_1  = R::pnorm(b, mu, sigma, true, false);
  double l1_2  = R::pnorm(a, mu, sigma, true, false);

  double l2_1  = R::pnorm(-a, mu, sigma, true, false);
  double l2_2  = R::pnorm(-b, mu, sigma, true, false);

  NumericVector L = log(l1 + l2) - log(l1_1-l1_2+l2_1-l2_2);
  return L;
}
// [[Rcpp::export(.tdist_lpdf)]]
NumericVector tdist_lpdf(NumericVector x, double mu, double df, double a, double b) {
  NumericVector l1    = Rcpp::dnt( x, df, mu, false);
  NumericVector l2    = Rcpp::dnt(-x, df, mu, false);
  
  double l1_1  = R::pnt(b, df, mu, true, false);
  double l1_2  = R::pnt(a, df, mu, true, false);
  
  double l2_1  = R::pnt(-a, df, mu, true, false);
  double l2_2  = R::pnt(-b, df, mu, true, false);
  
  NumericVector L = log(l1 + l2) - log(l1_1-l1_2+l2_1-l2_2);
  return L;
}
// [[Rcpp::export(.zdist_pdf)]]
NumericVector zdist_pdf(NumericVector x, double mu, double sigma, double a, double b) {
  NumericVector L = zdist_lpdf(x, mu, sigma, a, b);
  return exp(L);
}
// [[Rcpp::export(.zdist_cens_lpdf)]]
double zdist_cens_lpdf(double lb, double ub, double mu, double sigma, double a, double b) {
  
  double l1  = R::pnorm( ub, mu, sigma, true, false) - R::pnorm( lb, mu, sigma, true, false);
  double l2  = R::pnorm(-lb, mu, sigma, true, false) - R::pnorm(-ub, mu, sigma, true, false);

  double l1_1  = R::pnorm(b, mu, sigma, true, false);
  double l1_2  = R::pnorm(a, mu, sigma, true, false);
  
  double l2_1  = R::pnorm(-a, mu, sigma, true, false);
  double l2_2  = R::pnorm(-b, mu, sigma, true, false);
  
  double L = log(l1 + l2) - log(l1_1-l1_2+l2_1-l2_2);
  return L;
}
// [[Rcpp::export(.tdist_pdf)]]
NumericVector tdist_pdf(NumericVector x, double mu, double df, double a, double b) {
  NumericVector L = tdist_lpdf(x, mu, df, a, b);
  return exp(L);
}
NumericVector trunc_normal_lpdf(NumericVector x, double mu, double sigma, double a, double b) {
  double lccdf = R::pnorm(a, mu, sigma, false,true);
  double lcdf  = R::pnorm(b, mu, sigma, true, true);

  NumericVector L = dnorm(x, mu, sigma, true);

  L = L - lccdf - lcdf;
  return L;
}

double trunc_normal_E(double mu, double sigma, double a, double b) {

  double alfa = (a-mu)/sigma;
  double beta = (b-mu)/sigma;

  double E = mu + sigma*( (R::dnorm(alfa, 0, 1, false) - R::dnorm(beta, 0, 1, false)) / (R::pnorm(beta, 0, 1, true, false) - R::pnorm(alfa, 0, 1, true, false))  );

  return E;
}

// [[Rcpp::export(.dirichlet_rng)]]
NumericVector dirichlet_rng(NumericVector alpha){
  NumericVector random_value (alpha.length());

  for(int i = 0; i < alpha.length(); i++){
    random_value[i] = R::rgamma(alpha[i], 1);
  }
  random_value = random_value/sum(random_value);

  return random_value;
}
NumericVector random_mu(NumericVector mu_alpha, double mu_max){
  NumericVector new_mu  = (mu_alpha.size());
  NumericVector weights = cumsum(dirichlet_rng(mu_alpha));
  weights = 1 - weights;
  for(int i = 0; i < new_mu.length(); i++){
    new_mu[new_mu.length() - (i+1)] = weights[i];
  }
  new_mu = new_mu * mu_max;
  return new_mu;
}
NumericMatrix compute_u_log_lik(NumericVector x, NumericVector mu, NumericVector sigma, double a, double b){
  NumericMatrix ll(x.size(), mu.size());

  for(int k = 0; k < mu.size(); k++){
    ll(_,k) = zdist_lpdf(x,mu[k],sigma[k],a,b);
  }
  return ll;
}
NumericMatrix compute_u_log_lik_c(NumericVector x, NumericVector lb, NumericVector ub, NumericVector mu, NumericVector sigma, double a, double b){
  
  NumericMatrix ll_o(mu.size(), x.size());
  NumericMatrix ll_c(mu.size(), lb.size());
  
  for(int k = 0; k < mu.size(); k++){
    ll_o(k,_) = zdist_lpdf(x,mu[k],sigma[k],a,b);
  }
  
  for(int k = 0; k < mu.size(); k++){
    for(int i = 0; i < lb.size(); i++){
      ll_c(k,i) = zdist_cens_lpdf(lb[i],ub[i],mu[k],sigma[k],a,b);  
    }
  }
  
  NumericMatrix ll = transpose(cbind(ll_o, ll_c));
  
  return ll;
}
NumericMatrix weight_u_log_lik(NumericMatrix ull, NumericVector theta){
  NumericMatrix ll(ull.nrow(), ull.ncol());

  for(int k = 0; k < ull.ncol(); k++){
    ll(_,k) = ull(_,k) + log(theta[k]);
  }
  return ll;
}
NumericMatrix compute_log_lik(NumericVector x, NumericVector mu, NumericVector sigma, double a, double b, NumericVector theta){
  NumericMatrix ll(x.size(), mu.size());

  for(int k = 0; k < mu.size(); k++){
    ll(_,k) = trunc_normal_lpdf(x,mu[k],sigma[k],a,b) + log(theta[k]);
  }
  return ll;
}
double sum_finite(NumericMatrix ll){
  double s = 0;

  for(int k = 0; k < ll.ncol(); k++){
    LogicalVector inf_check = is_infinite(ll(_,k));
    NumericVector temp_ll = ll( _ ,k);
    temp_ll = temp_ll[!inf_check];
    s += sum(temp_ll);
  }
  return s;
}
NumericMatrix exp_matrix(NumericMatrix ll){
  NumericMatrix l (ll.nrow(), ll.ncol());
  for(int k = 0; k < ll.ncol(); k++){
    l(_,k) = exp(ll(_,k));
  }
  return(l);
}
NumericVector compute_l_row_sum(NumericMatrix l){
  NumericVector l_row_sum (l.nrow());
  for(int i = 0; i < l.nrow(); i++){
    l_row_sum[i] = sum(l(i,_));
  }
  return(l_row_sum);
}
NumericMatrix compute_p(NumericMatrix l, NumericVector l_row_sum){

  NumericMatrix p (l.nrow(), l.ncol());
  for(int i = 0; i < l.nrow(); i++){
    p(i,_) = l(i,_) / l_row_sum[i];
  }
  return(p);
}
NumericVector update_theta(NumericMatrix p){
  NumericVector new_theta (p.ncol());

  for(int k = 0; k < p.ncol(); k++){
    new_theta[k] = sum(p(_,k)) / p.nrow();
  }
  return(new_theta);
}
NumericVector bound_mu(NumericVector mu, double lower, double upper){
  NumericVector new_mu = ifelse(mu < lower, lower, mu);
  new_mu = ifelse(mu > upper, upper, mu);
  return new_mu;
}
NumericVector update_mu(NumericMatrix p, NumericVector x, NumericVector mu, NumericVector sigma, double a, double b){
  NumericVector new_mu = mu;

  for(int k = 1; k < p.ncol(); k ++){
    new_mu[k] = sum( p(_,k) * x  ) / sum( p(_,k) ) - trunc_normal_E(0, sigma[k], a - mu[k] , b - mu[k]);
  }
  new_mu = bound_mu(new_mu, 0, b + 2);
  return(new_mu);
}
NumericVector select_x(NumericVector x, double a, double b){
  LogicalVector x_true1 = x > a;
  LogicalVector x_true2 = x < b;
  NumericVector x_new  = x[x_true1 & x_true2];
  return x_new;
}
double get_prop_high(NumericVector x, double select_sig, double b){
  
  double a = R::pnorm(select_sig/2, 0, 1, false, false);
  
  LogicalVector x_sig_true = x > a;
  NumericVector x_sig      = x[x_sig_true];
  
  LogicalVector x_high_true = x > b;
  NumericVector x_high      = x[x_high_true];
  
  double prop_high = (1.0 * x_high.length()) / (1.0 * x_sig.length());
  return prop_high;
}
double get_prop_high_cens(NumericVector x, double select_sig, double b, int n_censored){
  
  double a = R::pnorm(select_sig/2, 0, 1, false, false);
  
  LogicalVector x_sig_true = x > a;
  NumericVector x_sig      = x[x_sig_true];
  
  LogicalVector x_high_true = x > b;
  NumericVector x_high      = x[x_high_true];
  
  double prop_high = (1.0 * x_high.length()) / (1.0 * (x_sig.length() + n_censored));
  return prop_high;
}


// [[Rcpp::export(.zcurve_EM_fit_RCpp)]]
List zcurve_EM_fit_RCpp(NumericVector x, int type, NumericVector mu, NumericVector sigma, NumericVector theta, double a, double b, double sig_level,
                 int max_iter, double criterion) {

  double prop_high = get_prop_high(x, sig_level, b);
  x = select_x(x,a,b);
  
  
  NumericMatrix log_lik (x.size(), mu.size());
  NumericMatrix lik (x.size(), mu.size());
  NumericVector l_row_sum (mu.size());
  NumericMatrix p (x.size(), mu.size());
  NumericVector Q (max_iter+1);

  int i= 0;
  Q[i] = 0;

  do{
    // E-step
    log_lik   = compute_log_lik(x, mu, sigma, a, b, theta);
    lik       = exp_matrix(log_lik);
    l_row_sum = compute_l_row_sum(lik);
    p         = compute_p(lik,l_row_sum);

    // M-step
    theta = update_theta(p);
    if(type == 2){
      mu = update_mu(p, x, mu, sigma, a, b);
    }

    Q[i+1] = sum(log(l_row_sum));
    ++i;

  } while ((fabs(Q[i]-Q[i-1]) >= criterion) && (i < max_iter));

  List ret;
  ret["iter"]      = i;
  ret["Q"]         = Q[i];
  ret["mu"]        = mu;
  ret["weights"]   = theta;
  ret["sigma"]     = sigma;
  ret["prop_high"] = prop_high;

  return ret;
}
// [[Rcpp::export(.zcurve_EM_fit_fast_RCpp)]]
List zcurve_EM_fit_fast_RCpp(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector theta, double a, double b, double sig_level,
                 int max_iter, double criterion) {

  double prop_high = get_prop_high(x, sig_level, b);
  x = select_x(x,a,b);
  
  NumericMatrix log_lik (x.size(), mu.size());
  NumericMatrix lik (x.size(), mu.size());
  NumericVector l_row_sum (mu.size());
  NumericMatrix p (x.size(), mu.size());
  NumericVector Q (max_iter+1);

  int i= 0;
  Q[i] = 0;


  NumericMatrix u_log_lik =  compute_u_log_lik(x, mu, sigma, a, b);
  do{
    // E-step
    log_lik   = weight_u_log_lik(u_log_lik, theta);
    lik       = exp_matrix(log_lik);
    l_row_sum = compute_l_row_sum(lik);
    p         = compute_p(lik,l_row_sum);

    // M-step
    theta = update_theta(p);

    Q[i+1] = sum(log(l_row_sum));
    ++i;

  } while ((fabs(Q[i]-Q[i-1]) >= criterion) && (i < max_iter));

  List ret;
  ret["iter"]      = i;
  ret["Q"]         = Q[i];
  ret["mu"]        = mu;
  ret["weights"]   = theta;
  ret["sigma"]     = sigma;
  ret["prop_high"] = prop_high;

  return ret;
}
// [[Rcpp::export(.zcurve_EMc_fit_fast_RCpp)]]
List zcurve_EMc_fit_fast_RCpp(NumericVector x, NumericVector lb, NumericVector ub,
                              NumericVector mu, NumericVector sigma, NumericVector theta, double a, double b, double sig_level,
                             int max_iter, double criterion) {
  
  double prop_high = get_prop_high_cens(x, sig_level, b, lb.size());
  x = select_x(x,a,b);
  
  NumericMatrix log_lik (x.size(), mu.size());
  NumericMatrix lik (x.size(), mu.size());
  NumericVector l_row_sum (mu.size());
  NumericMatrix p (x.size(), mu.size());
  NumericVector Q (max_iter+1);
  
  int i= 0;
  Q[i] = 0;
  
  
  NumericMatrix u_log_lik =  compute_u_log_lik_c(x, lb, ub, mu, sigma, a, b);
  do{
    // E-step
    log_lik   = weight_u_log_lik(u_log_lik, theta);
    lik       = exp_matrix(log_lik);
    l_row_sum = compute_l_row_sum(lik);
    p         = compute_p(lik,l_row_sum);
    
    // M-step
    theta = update_theta(p);
    
    Q[i+1] = sum(log(l_row_sum));
    ++i;
    
  } while ((fabs(Q[i]-Q[i-1]) >= criterion) && (i < max_iter));
  
  List ret;
  ret["iter"]      = i;
  ret["Q"]         = Q[i];
  ret["mu"]        = mu;
  ret["weights"]   = theta;
  ret["sigma"]     = sigma;
  ret["prop_high"] = prop_high;
  
  return ret;
}
// [[Rcpp::export(.zcurve_EM_start_RCpp)]]
List zcurve_EM_start_RCpp(NumericVector x, int type, int K,
                   NumericVector mu, NumericVector sigma,NumericVector mu_alpha, double mu_max,
                   NumericVector theta_alpha,
                   double a, double b, double sig_level,
                   int fit_reps, int max_iter, double criterion) {

  NumericMatrix mu_reps        (fit_reps, K);
  NumericMatrix weights_reps   (fit_reps, K);
  IntegerVector iter_reps      (fit_reps);
  NumericVector Q_reps         (fit_reps);
  NumericVector prop_high_reps (fit_reps);

  NumericVector temp_theta  (K);
  NumericVector temp_mu     (K);

  NumericVector new_mu      (K);
  NumericVector new_weights (K);
  int new_iter;
  double new_Q;
  double new_prop_high;

  for(int i = 0; i < fit_reps; i++){
    temp_theta = dirichlet_rng(theta_alpha);
    if(type == 1){
      temp_mu = mu;
    }else if(type == 2){
      temp_mu = random_mu(mu_alpha, mu_max);
    }

    List temp_fit = zcurve_EM_fit_RCpp(x, type, temp_mu, sigma, temp_theta, a, b, sig_level,
                                max_iter, criterion);

    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];

    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }

  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;

  return ret;
}

// [[Rcpp::export(.zcurve_EM_boot_RCpp)]]
List zcurve_EM_boot_RCpp(NumericVector x, int type,
                  NumericVector mu, NumericVector sigma, NumericVector theta,
                  double a, double b, double sig_level,
                  int bootstrap, int max_iter, double criterion){
  
  NumericMatrix mu_reps       (bootstrap, mu.size());
  NumericMatrix weights_reps  (bootstrap, mu.size());
  IntegerVector iter_reps     (bootstrap);
  NumericVector Q_reps        (bootstrap);
  NumericVector prop_high_reps (bootstrap);

  NumericVector temp_x;

  NumericVector new_mu      (mu.size());
  NumericVector new_weights (mu.size());
  int new_iter;
  double new_Q;
  double new_prop_high;

  for(int i = 0; i < bootstrap; i++){
    temp_x = sample(x, x.size(), true);

    List temp_fit = zcurve_EM_fit_RCpp(temp_x, type, mu, sigma, theta, a, b, sig_level,
                                max_iter, criterion);

    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];

    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }

  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;

  return ret;
}

// [[Rcpp::export(.zcurve_EM_start_fast_RCpp)]]
List zcurve_EM_start_fast_RCpp(NumericVector x, int K,
                   NumericVector mu, NumericVector sigma, NumericVector mu_alpha, double mu_max,
                   NumericVector theta_alpha,
                   double a, double b, double sig_level,
                   int fit_reps, int max_iter, double criterion) {

  NumericMatrix mu_reps        (fit_reps, K);
  NumericMatrix weights_reps   (fit_reps, K);
  IntegerVector iter_reps      (fit_reps);
  NumericVector Q_reps         (fit_reps);
  NumericVector prop_high_reps (fit_reps);

  NumericVector temp_theta  (K);
  NumericVector temp_mu     (K);

  NumericVector new_mu      (K);
  NumericVector new_weights (K);
  int new_iter;
  double new_Q;
  double new_prop_high;

  for(int i = 0; i < fit_reps; i++){
    temp_theta = dirichlet_rng(theta_alpha);
    temp_mu    = mu;

    List temp_fit = zcurve_EM_fit_fast_RCpp(x, temp_mu, sigma, temp_theta, a, b, sig_level,
                                max_iter, criterion);

    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];

    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }

  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;

  return ret;
}

// [[Rcpp::export(.zcurve_EM_boot_fast_RCpp)]]
List zcurve_EM_boot_fast_RCpp(NumericVector x,
                  NumericVector mu, NumericVector sigma, NumericVector theta,
                  double a, double b, double sig_level,
                  int bootstrap, int max_iter, double criterion){
  NumericMatrix mu_reps        (bootstrap, mu.size());
  NumericMatrix weights_reps   (bootstrap, mu.size());
  IntegerVector iter_reps      (bootstrap);
  NumericVector Q_reps         (bootstrap);
  NumericVector prop_high_reps (bootstrap);

  NumericVector temp_x;

  NumericVector new_mu      (mu.size());
  NumericVector new_weights (mu.size());
  int new_iter;
  double new_Q;
  double new_prop_high;

  for(int i = 0; i < bootstrap; i++){
    temp_x = sample(x, x.size(), true);

    List temp_fit = zcurve_EM_fit_fast_RCpp(temp_x, mu, sigma, theta, a, b, sig_level,
                                max_iter, criterion);

    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];

    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }

  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;

  return ret;
}

// [[Rcpp::export(.zcurve_EMc_start_fast_RCpp)]]
List zcurve_EMc_start_fast_RCpp(NumericVector x,  NumericVector lb, NumericVector ub, int K,
                               NumericVector mu, NumericVector sigma, NumericVector mu_alpha, double mu_max,
                               NumericVector theta_alpha,
                               double a, double b, double sig_level,
                               int fit_reps, int max_iter, double criterion) {
  
  NumericMatrix mu_reps        (fit_reps, K);
  NumericMatrix weights_reps   (fit_reps, K);
  IntegerVector iter_reps      (fit_reps);
  NumericVector Q_reps         (fit_reps);
  NumericVector prop_high_reps (fit_reps);
  
  NumericVector temp_theta  (K);
  NumericVector temp_mu     (K);
  
  NumericVector new_mu      (K);
  NumericVector new_weights (K);
  int new_iter;
  double new_Q;
  double new_prop_high;
  
  for(int i = 0; i < fit_reps; i++){
    temp_theta = dirichlet_rng(theta_alpha);
    temp_mu    = mu;
    
    List temp_fit = zcurve_EMc_fit_fast_RCpp(x, lb, ub, temp_mu, sigma, temp_theta, a, b, sig_level,
                                            max_iter, criterion);
    
    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];
    
    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }
  
  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;
  
  return ret;
}

// [[Rcpp::export(.zcurve_EMc_boot_fast_RCpp)]]
List zcurve_EMc_boot_fast_RCpp(NumericVector x, NumericVector lb, NumericVector ub,
                              NumericVector mu, NumericVector sigma, NumericVector theta,
                              double a, double b, double sig_level,
                              int bootstrap, int max_iter, double criterion){
  NumericMatrix mu_reps        (bootstrap, mu.size());
  NumericMatrix weights_reps   (bootstrap, mu.size());
  IntegerVector iter_reps      (bootstrap);
  NumericVector Q_reps         (bootstrap);
  NumericVector prop_high_reps (bootstrap);
  
  NumericVector temp_x;
  
  NumericVector new_mu      (mu.size());
  NumericVector new_weights (mu.size());
  int new_iter;
  double new_Q;
  double new_prop_high;
  
  for(int i = 0; i < bootstrap; i++){
    temp_x = sample(x, x.size(), true);
    
    List temp_fit = zcurve_EMc_fit_fast_RCpp(temp_x, lb, ub, mu, sigma, theta, a, b, sig_level,
                                            max_iter, criterion);
    
    new_mu        = temp_fit["mu"];
    new_weights   = temp_fit["weights"];
    new_iter      = temp_fit["iter"];
    new_Q         = temp_fit["Q"];
    new_prop_high = temp_fit["prop_high"];
    
    mu_reps(i,_)      = new_mu;
    weights_reps(i,_) = new_weights;
    iter_reps[i]      = new_iter;
    Q_reps[i]         = new_Q;
    prop_high_reps[i] = new_prop_high;
  }
  
  List ret;
  ret["iter"]      = iter_reps;
  ret["Q"]         = Q_reps;
  ret["mu"]        = mu_reps;
  ret["weights"]   = weights_reps;
  ret["prop_high"] = prop_high_reps;
  
  return ret;
}

/*
  double get_loglik(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector theta, double a, double b) {

  NumericMatrix log_lik   = compute_log_lik(x, mu, sigma, a, b, theta);
  NumericMatrix lik       = exp_matrix(log_lik);
  NumericVector l_row_sum = compute_l_row_sum(lik);
  NumericMatrix p         = compute_p(lik,l_row_sum);

  double  Q = sum(log(l_row_sum));

  return Q;
}
*/
