#include <cmath>         // for std::log or std::exp if used
#include <math.h>
#include "rng_utils.h"  // for safe_runif()
// Required headers
#include <RcppArmadillo.h>

using namespace Rcpp;


double p_inv_gamma(double dispersion,double shape,double rate){
  
  return(1- R::pgamma(1/dispersion,shape,1/rate,TRUE,FALSE));
}



double  q_inv_gamma(double p,double shape,double rate,double disp_upper,double disp_lower){
  double p_upp=p_inv_gamma(disp_upper,shape,rate);
  double p_low=p_inv_gamma(disp_lower,shape,rate);
  double p1=p_low+p*(p_upp-p_low);
  double p2=1-p1;
  return(1/ R::qgamma(p2,shape,1/rate,TRUE,FALSE));
}

namespace glmbayes {
namespace rng {

double rinvgamma_ct(double shape,double rate,double disp_upper,double disp_lower){
  double p= R::runif(0,1);
  return(q_inv_gamma(p,shape,rate,disp_upper,disp_lower));
}

}
}