// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Set_Grid
Rcpp::List Set_Grid(Rcpp::NumericMatrix GIndex, Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
RcppExport SEXP _glmbayes_Set_Grid(SEXP GIndexSEXP, SEXP cbarsSEXP, SEXP LintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type GIndex(GIndexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type cbars(cbarsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Lint(LintSEXP);
    rcpp_result_gen = Rcpp::wrap(Set_Grid(GIndex, cbars, Lint));
    return rcpp_result_gen;
END_RCPP
}
// setlogP
Rcpp::List setlogP(NumericMatrix logP, NumericVector NegLL, NumericMatrix cbars, NumericMatrix G3);
RcppExport SEXP _glmbayes_setlogP(SEXP logPSEXP, SEXP NegLLSEXP, SEXP cbarsSEXP, SEXP G3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type logP(logPSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NegLL(NegLLSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cbars(cbarsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G3(G3SEXP);
    rcpp_result_gen = Rcpp::wrap(setlogP(logP, NegLL, cbars, G3));
    return rcpp_result_gen;
END_RCPP
}
// f1_gamma
NumericVector f1_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_gamma(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_gamma
NumericVector f2_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f2_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_gamma(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f3_gamma
arma::mat f3_gamma(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f3_gamma(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_gamma(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f1_binomial_logit
NumericVector f1_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_binomial_logit(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_binomial_logit
NumericVector f2_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f2_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_binomial_logit(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f3_binomial_logit
arma::mat f3_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f3_binomial_logit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_binomial_logit(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f1_binomial_probit
NumericVector f1_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_binomial_probit(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_binomial_probit
NumericVector f2_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f2_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_binomial_probit(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f3_binomial_probit
arma::mat f3_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f3_binomial_probit(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_binomial_probit(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f1_binomial_cloglog
NumericVector f1_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_binomial_cloglog(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_binomial_cloglog
NumericVector f2_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f2_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f3_binomial_cloglog
arma::mat f3_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f3_binomial_cloglog(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f1_gaussian
NumericVector f1_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_gaussian(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_gaussian
NumericVector f2_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f2_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_gaussian(b, y, x, mu, P, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f3_gaussian
arma::mat f3_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f3_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_gaussian(b, y, x, mu, P, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f1_poisson
NumericVector f1_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f1_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_poisson(b, y, x, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_poisson
NumericVector f2_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f2_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f2_poisson(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// f3_poisson
arma::mat f3_poisson(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, int progbar);
RcppExport SEXP _glmbayes_f3_poisson(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(f3_poisson(b, y, x, mu, P, alpha, wt, progbar));
    return rcpp_result_gen;
END_RCPP
}
// glmbsim_cpp
Rcpp::List glmbsim_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_glmbsim_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< Function >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Envelope(EnvelopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(glmbsim_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// glmbenvelope_c
List glmbenvelope_c(NumericVector bStar, NumericMatrix A, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, std::string family, std::string link, int Gridtype, int n, bool sortgrid);
RcppExport SEXP _glmbayes_glmbenvelope_c(SEXP bStarSEXP, SEXP ASEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP, SEXP nSEXP, SEXP sortgridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type bStar(bStarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type Gridtype(GridtypeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type sortgrid(sortgridSEXP);
    rcpp_result_gen = Rcpp::wrap(glmbenvelope_c(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, sortgrid));
    return rcpp_result_gen;
END_RCPP
}
// glmbsim_NGauss_cpp
Rcpp::List glmbsim_NGauss_cpp(int n, NumericVector y, NumericMatrix x, NumericVector mu, NumericMatrix P, NumericVector offset2, NumericVector wt, double dispersion, Rcpp::List famfunc, Function f1, Function f2, Function f3, NumericVector start, std::string family, std::string link, int Gridtype);
RcppExport SEXP _glmbayes_glmbsim_NGauss_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP offset2SEXP, SEXP wtSEXP, SEXP dispersionSEXP, SEXP famfuncSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP startSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset2(offset2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type famfunc(famfuncSEXP);
    Rcpp::traits::input_parameter< Function >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< Function >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< Function >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type Gridtype(GridtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(glmbsim_NGauss_cpp(n, y, x, mu, P, offset2, wt, dispersion, famfunc, f1, f2, f3, start, family, link, Gridtype));
    return rcpp_result_gen;
END_RCPP
}
// glmbsim_Gauss_cpp
Rcpp::List glmbsim_Gauss_cpp(int n, NumericVector y, NumericMatrix x, NumericVector mu, NumericMatrix P, NumericVector offset2, NumericVector wt, double dispersion, Rcpp::List famfunc, Function f1, Function f2, Function f3, NumericVector start, std::string family, std::string link, int Gridtype);
RcppExport SEXP _glmbayes_glmbsim_Gauss_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP offset2SEXP, SEXP wtSEXP, SEXP dispersionSEXP, SEXP famfuncSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP startSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset2(offset2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type famfunc(famfuncSEXP);
    Rcpp::traits::input_parameter< Function >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< Function >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< Function >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type Gridtype(GridtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(glmbsim_Gauss_cpp(n, y, x, mu, P, offset2, wt, dispersion, famfunc, f1, f2, f3, start, family, link, Gridtype));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_glmbayes_Set_Grid", (DL_FUNC) &_glmbayes_Set_Grid, 3},
    {"_glmbayes_setlogP", (DL_FUNC) &_glmbayes_setlogP, 4},
    {"_glmbayes_f1_gamma", (DL_FUNC) &_glmbayes_f1_gamma, 5},
    {"_glmbayes_f2_gamma", (DL_FUNC) &_glmbayes_f2_gamma, 8},
    {"_glmbayes_f3_gamma", (DL_FUNC) &_glmbayes_f3_gamma, 8},
    {"_glmbayes_f1_binomial_logit", (DL_FUNC) &_glmbayes_f1_binomial_logit, 5},
    {"_glmbayes_f2_binomial_logit", (DL_FUNC) &_glmbayes_f2_binomial_logit, 8},
    {"_glmbayes_f3_binomial_logit", (DL_FUNC) &_glmbayes_f3_binomial_logit, 8},
    {"_glmbayes_f1_binomial_probit", (DL_FUNC) &_glmbayes_f1_binomial_probit, 5},
    {"_glmbayes_f2_binomial_probit", (DL_FUNC) &_glmbayes_f2_binomial_probit, 8},
    {"_glmbayes_f3_binomial_probit", (DL_FUNC) &_glmbayes_f3_binomial_probit, 8},
    {"_glmbayes_f1_binomial_cloglog", (DL_FUNC) &_glmbayes_f1_binomial_cloglog, 5},
    {"_glmbayes_f2_binomial_cloglog", (DL_FUNC) &_glmbayes_f2_binomial_cloglog, 8},
    {"_glmbayes_f3_binomial_cloglog", (DL_FUNC) &_glmbayes_f3_binomial_cloglog, 8},
    {"_glmbayes_f1_gaussian", (DL_FUNC) &_glmbayes_f1_gaussian, 5},
    {"_glmbayes_f2_gaussian", (DL_FUNC) &_glmbayes_f2_gaussian, 7},
    {"_glmbayes_f3_gaussian", (DL_FUNC) &_glmbayes_f3_gaussian, 7},
    {"_glmbayes_f1_poisson", (DL_FUNC) &_glmbayes_f1_poisson, 5},
    {"_glmbayes_f2_poisson", (DL_FUNC) &_glmbayes_f2_poisson, 8},
    {"_glmbayes_f3_poisson", (DL_FUNC) &_glmbayes_f3_poisson, 8},
    {"_glmbayes_glmbsim_cpp", (DL_FUNC) &_glmbayes_glmbsim_cpp, 12},
    {"_glmbayes_glmbenvelope_c", (DL_FUNC) &_glmbayes_glmbenvelope_c, 13},
    {"_glmbayes_glmbsim_NGauss_cpp", (DL_FUNC) &_glmbayes_glmbsim_NGauss_cpp, 16},
    {"_glmbayes_glmbsim_Gauss_cpp", (DL_FUNC) &_glmbayes_glmbsim_Gauss_cpp, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_glmbayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
