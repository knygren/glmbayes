// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Set_Grid
Rcpp::List Set_Grid(Rcpp::NumericMatrix GIndex, Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
RcppExport SEXP glmbayes_Set_Grid(SEXP GIndexSEXP, SEXP cbarsSEXP, SEXP LintSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type GIndex(GIndexSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type cbars(cbarsSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Lint(LintSEXP );
        Rcpp::List __result = Set_Grid(GIndex, cbars, Lint);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// setlogP
Rcpp::List setlogP(NumericMatrix logP, NumericVector NegLL, NumericMatrix cbars, NumericMatrix G3);
RcppExport SEXP glmbayes_setlogP(SEXP logPSEXP, SEXP NegLLSEXP, SEXP cbarsSEXP, SEXP G3SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type logP(logPSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type NegLL(NegLLSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type cbars(cbarsSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G3(G3SEXP );
        Rcpp::List __result = setlogP(logP, NegLL, cbars, G3);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_binomial_logit
NumericVector f1_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_binomial_logit(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_binomial_logit
NumericVector f2_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_binomial_logit(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_binomial_logit
arma::mat f3_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_binomial_logit(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_binomial_probit
NumericVector f1_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_binomial_probit(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_binomial_probit
NumericVector f2_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_binomial_probit(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_binomial_probit
arma::mat f3_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_binomial_probit(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_binomial_cloglog
NumericVector f1_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_binomial_cloglog(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_binomial_cloglog
NumericVector f2_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_binomial_cloglog(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_binomial_cloglog
arma::mat f3_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_binomial_cloglog(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_gamma
NumericVector f1_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_gamma(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_gamma
NumericVector f2_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_gamma(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_gamma
arma::mat f3_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_gamma(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_gaussian
NumericVector f1_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_gaussian(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_gaussian
NumericVector f2_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_gaussian(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_gaussian
arma::mat f3_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_gaussian(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f1_poisson
NumericVector f1_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP glmbayes_f1_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f1_poisson(b, y, x, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f2_poisson
NumericVector f2_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f2_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        NumericVector __result = f2_poisson(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// f3_poisson
arma::mat f3_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP glmbayes_f3_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        arma::mat __result = f3_poisson(b, y, x, mu, P, alpha, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// glmbsim_cpp
Rcpp::List glmbsim_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::CharacterVector family, Rcpp::CharacterVector link);
RcppExport SEXP glmbayes_glmbsim_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP familySEXP, SEXP linkSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        Rcpp::traits::input_parameter< Function >::type f2(f2SEXP );
        Rcpp::traits::input_parameter< Rcpp::List >::type Envelope(EnvelopeSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP );
        Rcpp::List __result = glmbsim_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
