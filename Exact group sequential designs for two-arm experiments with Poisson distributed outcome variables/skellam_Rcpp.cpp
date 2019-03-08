#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector dskellam_rcpp(IntegerVector x, double lambda1, double lambda2) {
  int           len_x = x.length();
  NumericVector f(len_x);
  if ((lambda1 == lambda2) && (lambda1 != 0)) {
    double K2,
           s,
           spd,
           sum_lambda          = lambda1 + lambda2,
           log_factor_2        = log(1/sum_lambda),
           add_factor_2        = pow(sum_lambda, 2),
           pi_factor           = sqrt(2*M_PI);
    if (x[0] >= 0) {
      //Rcpp::Rcout << "...probs in pos_t2" << std::endl;
      for (int i = 0; i < len_x; i++) {
        //Rcpp::Rcout << "...probs in f" << std::endl;
        f[i]                   = R::bessel_i(sum_lambda, x[i], 2);
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*(exp(s) + exp(-s));
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*(exp(s) - exp(-s)))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
    else if (x[len_x - 1] <= 0) {
      //Rcpp::Rcout << "...probs in neg_t2" << std::endl;
      NumericVector abs_x      = abs(x);
      for (int i = 0; i < len_x; i++) {
        //Rcpp::Rcout << "...probs in f" << std::endl;
        f[i]                   = R::bessel_i(sum_lambda, abs_x[i], 2);
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*(exp(s) + exp(-s));
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*(exp(s) - exp(-s)))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
    else {
      NumericVector abs_x      = abs(x),
                    bessels(max(abs_x) + 1);
      for (int i = 0; i < len_x; i++) {
       // Rcpp::Rcout << "...probs in posneg" << std::endl;
        if (bessels[abs_x[i]] == 0) {
          bessels[abs_x[i]]    = R::bessel_i(sum_lambda, abs_x[i], 2);
        }
        //Rcpp::Rcout << "...probs in f" << std::endl;
        f[i]                   = bessels[abs_x[i]];
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*(exp(s) + exp(-s));
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*(exp(s) - exp(-s)))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
  }
  else if ((lambda1 != 0) && (lambda2 != 0)) {
    double        K2,
                  s,
                  spd,
                  sum_lambda   = lambda1 + lambda2,
                  log_factor   = 0.5*log(lambda1/lambda2),
                  log_factor_2 = log(0.5/lambda2),
                  sqrt_factor  = 2*sqrt(lambda2*lambda1),
                  add_factor   = sqrt_factor - sum_lambda,
                  add_factor_2 = pow(sqrt_factor, 2),
                  pi_factor    = sqrt(2*M_PI);
    if (x[0] >= 0) {
      //Rcpp::Rcout << "...probs in pos" << std::endl;
      for (int i = 0; i < len_x; i++) {
        //Rcpp::Rcout << "...probs in f" << std::endl;
        f[i]                   = R::bessel_i(sqrt_factor, x[i], 2)*
                                   exp(add_factor + x[i]*log_factor);
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*exp(s) + lambda2*exp(-s);
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*exp(s) - lambda2*exp(-s))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
    else if (x[len_x - 1] <= 0) {
      //Rcpp::Rcout << "...probs in neg" << std::endl;
      NumericVector abs_x      = abs(x);
      for (int i = 0; i < len_x; i++) {
        //Rcpp::Rcout << "...probs in f" << std::endl;
        f[i]                   = R::bessel_i(sqrt_factor, abs_x[i], 2)*
                                   exp(add_factor + x[i]*log_factor);
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*exp(s) + lambda2*exp(-s);
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*exp(s) - lambda2*exp(-s))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
    else {
      //Rcpp::Rcout << "...probs in mix" << std::endl;
      NumericVector abs_x      = abs(x),
                    bessels(max(abs_x) + 1);
      for (int i = 0; i < len_x; i++) {
        //Rcpp::Rcout << "...probs in f" << std::endl;
        if (bessels[abs_x[i]] == 0) {
          bessels[abs_x[i]]    = R::bessel_i(sqrt_factor, abs_x[i], 2);
        }
        f[i]                   = bessels[abs_x[i]]*
                                   exp(add_factor + x[i]*log_factor);
        if (f[i] < 1e-308) {
          //Rcpp::Rcout << "...probs in small f" << std::endl;
          s                    = log_factor_2 + log(x[i] + sqrt(pow(x[i], 2) +
                                                      add_factor_2));
          K2                   = lambda1*exp(s) + lambda2*exp(-s);
          spd                  = exp(K2 - sum_lambda - x[i]*s)/
                                   (pi_factor*sqrt(K2));
          if ((spd > 1e-308) && (spd != R_PosInf)) {
            //Rcpp::Rcout << "...probs in spd" << std::endl;
            f[i]               =
              spd*(2 + 0.125*(1 - 5*pow((lambda1*exp(s) - lambda2*exp(-s))/K2,
                                        2)/3)/K2)*0.5;
          }
        }
      }
    }
  }
  else if ((lambda1 == 0) && (lambda2 == 0)) {
    for (int i = 0; i < len_x; i++) {
      if (x[i] == 0) {
        f[i]                   = 1;
      }
    }
  }
  else if (lambda2 == 0) {
    for (int i = 0; i < len_x; i++) {
      f[i]                     = R::dpois(x[i], lambda1, 0);
    }
  }
  else {
    for (int i = 0; i < len_x; i++) {
      f[i]                     = R::dpois(-x[i], lambda2, 0);
    }
  }
  return f;
}

// [[Rcpp::export]]
NumericVector pskellam_rcpp(IntegerVector x, double lambda1, double lambda2,
                            int lower_tail) {
  int           len_x       = x.length();
  double        two_lambda1 = 2*lambda1,
                two_lambda2 = 2*lambda2;
  NumericVector cf(len_x);
   if (lower_tail == 1) {
     //Rcpp::Rcout << "...probs in lt" << std::endl;
     for (int i = 0; i < len_x; i++) {
       if (x[i] < 0) {
         //Rcpp::Rcout << "...probs in xi<0" << std::endl;
         cf[i]              = R::pnchisq(two_lambda2, -2*x[i], two_lambda1,
                                           1, 0);
       }
       else {
         //Rcpp::Rcout << "...probs in xi>=0" << std::endl;
         cf[i]              = R::pnchisq(two_lambda1, 2*x[i] + 2, two_lambda2,
                                         0, 0);
         //Rcpp::Rcout << "...probs in xi>=0:done" << std::endl;
       }
     }
   }
   else {
     //Rcpp::Rcout << "...probs in nlt" << std::endl;
     for (int i = 0; i < len_x; i++) {
       if (x[i] < 0) {
         //Rcpp::Rcout << "...probs in xi<0" << std::endl;
         cf[i]              = R::pnchisq(two_lambda2, -2*x[i], two_lambda1, 0,
                                         0);
       }
       else {
         //Rcpp::Rcout << "...probs in xi>=0" << std::endl;
         cf[i]              = R::pnchisq(two_lambda1, 2*x[i] + 2, two_lambda2,
                                         1, 0);
         //Rcpp::Rcout << "...probs in xi>=0:done" << std::endl;
       }
     }
   }
   //Rcpp::Rcout << "...probs in adjustmentpart" << std::endl;
   double xm,
          s,
          K2,
          lambda1_mod,
          lambda2_mod,
          u2,
          w2,
          xe;
   IntegerVector x_mod;
   if (lower_tail == 0) {
     x_mod                  = -x - 1;
     lambda1_mod            = lambda2;
     lambda2_mod            = lambda1;
   }
   else {
     IntegerVector x_mod    = x;
     lambda1_mod     = lambda1;
     lambda2_mod     = lambda2;
   }
   double g16               = (lambda1_mod - lambda2_mod)/
                                (6*pow(lambda1_mod + lambda2_mod, 1.5)),
     lambda_factor          = 4*lambda1_mod*lambda2_mod,
     log_factor             = log(0.5/lambda2_mod),
     sum_lambda             = lambda1_mod + lambda2_mod,
     sqrt_sum_lambda        = sqrt(sum_lambda),
     lambda_factor_2        = (lambda1_mod - lambda2_mod)/sqrt(sum_lambda);
   for (int i = 0; i < len_x; i++) {
     if (cf[i] < 1e-308) {
       xm                   = -x[i] - 0.5;
       s                    = log_factor +
                                log(xm + sqrt(pow(xm, 2) + lambda_factor));
       K2                   = lambda2_mod*exp(s) + lambda1_mod*exp(-s);
       u2                   = 2*sinh(0.5*s)*sqrt(K2);
       w2                   = R::sign(s)*sqrt(2*(s*xm - K2 + sum_lambda));
       xe                   = xm/sqrt_sum_lambda + lambda_factor_2;
       if (abs(xe) < 1e-4) {
         cf[i]              = R::pnorm(-xe, 0, 1, 1, 0) +
                                R::dnorm(xe, 0, 1, 0)*g16*(1 - pow(xe, 2));
       }
       else {
         cf[i]              = R::pnorm(-w2, 0, 1, 1, 0) -
                                R::dnorm(w2, 0, 1, 0)*(1/w2 - 1/u2);
       }
     }
   }
   return cf;
}

// [[Rcpp::export]]
double power_exact_rcpp(double lambda1, double lambda2, int K, int n,
                        IntegerVector a, IntegerVector r) {
  double P                  = pskellam_rcpp(IntegerVector::create(r[0] - 1),
                                            n*lambda1, n*lambda2, 0)[0];
  if (K > 1) {
    NumericVector pmf_k_old = dskellam_rcpp(seq(a[0] + 1, r[0] - 1), n*lambda1,
                                            n*lambda2);
    if (K > 2) {
      for (int k = 2; k <= K - 1; k++) {
        NumericVector pmf_k_new(r[k - 1] - a[k - 1] - 1);
        for (int l = a[k - 1] + 1; l <= r[k - 1] - 1; l++) {
          pmf_k_new         = pmf_k_new +
                                pmf_k_old[l - a[k - 1] - 1]*
                                  dskellam_rcpp(seq(a[k - 1] + 1,
                                                    r[k - 1] - 1) - l,
                                                n*lambda1, n*lambda2);
          P                += pmf_k_old[l - a[k - 1] - 1]*
                                pskellam_rcpp(r[k - 1] - l - 1, n*lambda1,
                                              n*lambda2, 0)[0];
        }
        pmf_k_old           = pmf_k_new;
      }
    }
    for (int l = a[K - 1] + 1; l <= r[K - 1] - 1; l++) {
      P                    += pmf_k_old[l - a[K - 1] - 1]*
                                pskellam_rcpp(r[K] - l - 1, n*lambda1,
                                              n*lambda2, 0)[0];
    }
  }
  return P;
}

// [[Rcpp::export]]
NumericVector exact_max_typeI(int K, double alpha, int n, IntegerVector a,
                              IntegerVector r, NumericVector lambda_null) {
  int           iter    = 0;
  double        f_v,
  golden  = 0.5*(3 - pow(5, 0.5)),
  x_left  = lambda_null[0],
  x_right = lambda_null[1],
  v       = x_left + golden*(x_right - x_left);
  f_v                 = -power_exact_rcpp(v, v, K, n, a, r);
  NumericVector output(4);
  double f_z,
  z              = 0.5*sum(lambda_null);
  f_z                 = -power_exact_rcpp(z, z, K, n, a, r);
  double f_u,
  midpoint,
  u,
  w              = v,
  f_w            = f_v,
  w_left         = z - x_left,
  w_right        = x_right - z,
  tol            = 1e-4*z,
  two_tol        = 2*tol,
  d              = 0,
  ed             = 0,
  p              = 0,
  q              = 0,
  rb              = 0;
  while (iter < 100) {
    midpoint            = 0.5*(x_left + x_right);
    if (abs(z - midpoint) <= two_tol - 0.5*(x_right - x_left)) {
      output            = NumericVector::create(z, -f_z, 0, iter);
      return output;
    }
    if (abs(ed) > tol) {
      rb                 = (z - w)*(f_z - f_v);
      q                 = (z - v)*(f_z - f_w);
      p                 = (z - v)*q - (z - w)*rb;
      q                 = 2*(q - rb);
      if (q > 0) {
        p               = -p;
      }
      else {
        q               = -q;
      }
      rb                 = ed;
      ed                = d;
    }
    if ((abs(p) < abs(0.5*q*rb)) && (p < q*w_left) && (p < q*w_right)) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) || (x_right - u < two_tol)) {
        d               = (z < midpoint ? tol : -tol);
      }
    }
    else {
      ed                = (z < midpoint ? x_right - z : -(z - x_left));
      d                 = golden*ed;
    }
    if (abs(d) >= tol) {
      u                 = z + d;
    }
    else {
      u                 = z + (d > 0 ? tol : -tol);
    }
    f_u                 = -power_exact_rcpp(u, u, K, n, a, r);
    if (f_u <= f_z) {
      if (u < z) {
        x_right         = z;
      }
      else {
        x_left          = z;
      }
      v                 = w;
      f_v               = f_w;
      w                 = z;
      f_w               = f_z;
      z                 = u;
      f_z               = f_u;
    }
    else {
      if (u < z) {
        x_left          = u;
      }
      else {
        x_right         = u;
      }
      if ((f_u <= f_w) || (w == z)) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if ((f_u <= f_v) || (v == z) || (v == w)) {
        v               = u;
        f_v             = f_u;
      }
    }
    w_left              = z - x_left;
    w_right             = x_right - z;
    iter++;
  }
  output                = NumericVector::create(z, -f_z, 1, iter);
  return output;
}

// [[Rcpp::export]]
NumericVector exact_min_power(int K, double beta, double delta, int n,
                              IntegerVector a, IntegerVector r,
                              NumericVector lambda_alt) {
  int           iter    = 0;
  double        f_v,
  golden  = 0.5*(3 - pow(5, 0.5)),
  x_left  = lambda_alt[0],
                       x_right = lambda_alt[1],
                                            v       = x_left + golden*(x_right - x_left);
  f_v                 = power_exact_rcpp(v, v - delta, K, n, a, r);
  NumericVector output(4);
  double f_z,
  z              = 0.5*sum(lambda_alt);
  f_z                 = power_exact_rcpp(z, z - delta, K, n, a, r);
  double f_u,
  midpoint,
  u,
  w              = v,
  f_w            = f_v,
  w_left         = z - x_left,
  w_right        = x_right - z,
  tol            = 1e-4*z,
  two_tol        = 2*tol,
  d              = 0,
  ed             = 0,
  p              = 0,
  q              = 0,
  rb              = 0;
  while (iter < 100) {
    midpoint            = 0.5*(x_left + x_right);
    if (abs(z - midpoint) <= two_tol - 0.5*(x_right - x_left)) {
      output            = NumericVector::create(z, f_z, 0, iter);
      return output;
    }
    if (abs(ed) > tol) {
      rb                 = (z - w)*(f_z - f_v);
      q                 = (z - v)*(f_z - f_w);
      p                 = (z - v)*q - (z - w)*rb;
      q                 = 2*(q - rb);
      if (q > 0) {
        p               = -p;
      }
      else {
        q               = -q;
      }
      rb                 = ed;
      ed                = d;
    }
    if ((abs(p) < abs(0.5*q*rb)) && (p < q*w_left) && (p < q*w_right)) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) || (x_right - u < two_tol)) {
        d               = (z < midpoint ? tol : -tol);
      }
    }
    else {
      ed                = (z < midpoint ? x_right - z : -(z - x_left));
      d                 = golden*ed;
    }
    if (abs(d) >= tol) {
      u                 = z + d;
    }
    else {
      u                 = z + (d > 0 ? tol : -tol);
    }
    f_u                 = power_exact_rcpp(u, u - delta, K, n, a, r);
    if (f_u <= f_z) {
      if (u < z) {
        x_right         = z;
      }
      else {
        x_left          = z;
      }
      v                 = w;
      f_v               = f_w;
      w                 = z;
      f_w               = f_z;
      z                 = u;
      f_z               = f_u;
    }
    else {
      if (u < z) {
        x_left          = u;
      }
      else {
        x_right         = u;
      }
      if ((f_u <= f_w) || (w == z)) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if ((f_u <= f_v) || (v == z) || (v == w)) {
        v               = u;
        f_v             = f_u;
      }
    }
    w_left              = z - x_left;
    w_right             = x_right - z;
    iter++;
  }
  output                = NumericVector::create(z, f_z, 1, iter);
  return output;
}

// [[Rcpp::export]]
NumericMatrix exact_des_two_stage_cpp(double alpha, double beta,
                                      double delta, int min_n, int max_n,
                                      double min_lambda_null,
                                      double max_lambda_null,
                                      double min_lambda_alt,
                                      double max_lambda_alt,
                                      double lambda_ess, int summary) {
  int           a_n,
                r_n,
                counter      = 0;
  double        ess0,
                ess1,
                typeI_min,
                typeI_max,
                typeI_check_min,
                typeI_check_max,
                typeI1_min,
                typeII_min,
                typeI1_max,
                typeII_max,
                typeII_check_min,
                typeII_check_max,
                typeII1_min,
                typeII1_max;
  IntegerVector poss_ar,
                poss_dp;
  NumericVector dskellam_ess,
                dskellam_power_min,
                dskellam_power_max,
                dskellam_typeI_min,
                dskellam_typeI_max,
                pskellam_a_n_min,
                pskellam_a_n_2_min,
                pskellam_a_n_max,
                pskellam_a_n_2_max,
                pskellam_ess0_a,
                pskellam_ess0_r,
                pskellam_ess1_a,
                pskellam_ess1_r,
                pskellam_r_n_min,
                pskellam_r_n_2_min,
                pskellam_r_n_max,
                pskellam_r_n_2_max;
  NumericMatrix feasible_designs(20000000, 10);
  for (int n = min_n; n <= max_n; n++) {
    if (summary == 1) {
      Rcpp::Rcout << "...currently analysing designs with n = " << n << "..." <<
        std::endl;
    }
    //Rcpp::Rcout << "...probs in a_n" << std::endl;
    a_n                      = 1;
    typeII_check_max         = 1,
    typeII_check_min         = 1;
    while ((typeII_check_max > 1e-7) || (typeII_check_min > 1e-7)) {
      a_n--;
      typeII_check_max           = pskellam_rcpp(IntegerVector::create(a_n),
                                             n*max_lambda_alt,
                                             n*(max_lambda_alt - delta), 1)[0];
      typeII_check_min           = pskellam_rcpp(IntegerVector::create(a_n),
                                             n*min_lambda_alt,
                                             n*(min_lambda_alt - delta), 1)[0];
    }
    //Rcpp::Rcout << "...probs in r_n" << std::endl;
    r_n                      = -1;
    typeI_check_max              = 1;
    typeI_check_min         = 1;
    while ((typeI_check_max > 1e-7) || (typeI_check_min > 1e-7)) {
      r_n++;
      typeI_check_max            = pskellam_rcpp(IntegerVector::create(r_n),
                                             n*max_lambda_null, n*max_lambda_null,
                                             0)[0];
      typeI_check_min            = pskellam_rcpp(IntegerVector::create(r_n),
                                                 n*min_lambda_null, n*min_lambda_null,
                                                 0)[0];
    }
    r_n++;
    //Rcpp::Rcout << "...probs in a_n r_n" << std::endl;
    if (a_n < r_n - 1) {
      //Rcpp::Rcout << "..." << a_n << "..." << r_n << std::endl;
      poss_ar                = seq(a_n, r_n);
      poss_dp                = seq(2*a_n - r_n, 2*r_n - a_n);
      //Rcpp::Rcout << "...probs in here2" << std::endl;
      pskellam_a_n_max           = pskellam_rcpp(poss_ar, n*max_lambda_alt,
                                             n*(max_lambda_alt - delta), 1);
      pskellam_a_n_min           = pskellam_rcpp(poss_ar, n*min_lambda_alt,
                                                 n*(min_lambda_alt - delta), 1);
      //Rcpp::Rcout << "...probs in here3" << std::endl;
      pskellam_r_n_max           = pskellam_rcpp(poss_ar - 1, n*max_lambda_null,
                                             n*max_lambda_null, 0);
      pskellam_r_n_min           = pskellam_rcpp(poss_ar - 1, n*min_lambda_null,
                                                 n*min_lambda_null, 0);
      //Rcpp::Rcout << "...probs in here4" << std::endl;
      pskellam_a_n_2_max         = pskellam_rcpp(poss_dp, n*max_lambda_alt,
                                             n*(max_lambda_alt - delta), 1);
      pskellam_a_n_2_min         = pskellam_rcpp(poss_dp, n*min_lambda_alt,
                                                 n*(min_lambda_alt - delta), 1);
      //Rcpp::Rcout << "...probs in here5" << std::endl;
      pskellam_r_n_2_max         = pskellam_rcpp(poss_dp - 1, n*max_lambda_null,
                                             n*max_lambda_null, 0);
      pskellam_r_n_2_min         = pskellam_rcpp(poss_dp - 1, n*min_lambda_null,
                                             n*min_lambda_null, 0);
      //Rcpp::Rcout << "...probs in here6" << std::endl;
      pskellam_ess0_a        = pskellam_rcpp(poss_ar, n*lambda_ess,
                                             n*lambda_ess, 1);
      //Rcpp::Rcout << "...probs in here7" << std::endl;
      pskellam_ess0_r        = pskellam_rcpp(poss_ar - 1, n*lambda_ess,
                                             n*lambda_ess, 0);
      //Rcpp::Rcout << "...probs in here8" << std::endl;
      pskellam_ess1_a        = pskellam_rcpp(poss_ar, n*lambda_ess,
                                             n*(lambda_ess - delta), 1);
      //Rcpp::Rcout << "...probs in here9" << std::endl;
      pskellam_ess1_r        = pskellam_rcpp(poss_ar - 1, n*lambda_ess,
                                             n*(lambda_ess - delta), 0);
      //Rcpp::Rcout << "...probs in here10" << std::endl;
      dskellam_typeI_max         = dskellam_rcpp(poss_ar, n*max_lambda_null,
                                             n*max_lambda_null);
      dskellam_typeI_min         = dskellam_rcpp(poss_ar, n*min_lambda_null,
                                                 n*min_lambda_null);
      //Rcpp::Rcout << "...probs in here11" << std::endl;
      dskellam_power_max         = dskellam_rcpp(poss_ar, n*max_lambda_alt,
                                             n*(max_lambda_alt - delta));
      dskellam_power_min         = dskellam_rcpp(poss_ar, n*min_lambda_alt,
                                             n*(min_lambda_alt - delta));
      //Rcpp::Rcout << "...probs in here12" << std::endl;
      for (int a1 = a_n; a1 < r_n - 1; a1++) {
        typeII1_max              = pskellam_a_n_max[a1 - a_n];
        typeII1_min              = pskellam_a_n_min[a1 - a_n];
        if ((typeII1_max < beta) && (typeII1_min < beta)) {
          for (int r1 = a1 + 2; r1 <= r_n; r1++) {
            typeI1_max           = pskellam_r_n_max[r1 - a_n];
            typeI1_min           = pskellam_r_n_min[r1 - a_n];
            ess0             = 2*n*(2 - pskellam_ess0_a[a1 - a_n] -
                                      pskellam_ess0_r[r1 - a_n]),
            ess1             = 2*n*(2 - pskellam_ess1_a[a1 - a_n] -
                                      pskellam_ess1_r[r1 - a_n]);
            if ((typeI1_max < alpha) && (typeI1_min < alpha)) {
              for (int r2 = a1 + a_n; r2 <= r1 + r_n; r2++) {
                //if (n == 98 && a1 == 7 && r1 == 94 && r2 == 131) {
                //  Rcpp::Rcout << "...ess0 = " << ess0 << "...ess1 = " << ess1 << std::endl;
                //}
                typeI_max       = typeI1_max;
                typeII_max      = typeII1_max;
                typeI_min       = typeI1_min;
                typeII_min      = typeII1_min;
                for (int l = a1 + 1; l <= r1 - 1; l++) {
                 typeI_max    += dskellam_typeI_max[l - a_n]*
                                 pskellam_r_n_2_max[r2 - l - (2*a_n - r_n)];
                  typeII_max   += dskellam_power_max[l - a_n]*
                                 pskellam_a_n_2_max[r2 - l - 1 - (2*a_n - r_n)];
                  typeI_min    += dskellam_typeI_min[l - a_n]*
                    pskellam_r_n_2_min[r2 - l - (2*a_n - r_n)];
                  typeII_min   += dskellam_power_min[l - a_n]*
                    pskellam_a_n_2_min[r2 - l - 1 - (2*a_n - r_n)];
                  if ((typeI_max > alpha) || (typeII_max > beta) || (typeI_min > alpha) || (typeII_min > beta)) {
                    break;
                  }
                }
                //if (n == 98 && a1 == 7 && r1 == 94 && r2 == 131) {
                //  Rcpp::Rcout << "...typeI = " << typeI << "...typeII = " << typeII << std::endl;
                //}
                if ((typeI_max <= alpha) && (typeII_max < beta) && (typeI_min <= alpha) && (typeII_min < beta)) {
                  //NumericVector maxtypeI = exact_max_typeI(2,  0.05, n,
                  //                                         IntegerVector::create(a1, r2 - 1),
                  //                                         IntegerVector::create(r1, r2),
                  //                                         NumericVector::create(0.25*max_lambda_null, max_lambda_null));
                  feasible_designs(counter, _) =
                    NumericVector::create(n, a1, r1, r2, typeI_min, typeI_max, 1 - typeII_min, 1 - typeII_max, ess0, ess1);
                  counter++;
                }
                else if ((typeII_min >= beta) || (typeII_max >= beta)) {
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output                    =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 9));
  return output;
}
