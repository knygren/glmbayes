// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EnvelopeBuild_c
List EnvelopeBuild_c(NumericVector bStar, NumericMatrix A, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, std::string family, std::string link, int Gridtype, int n, bool sortgrid);
RcppExport SEXP _glmbayes_EnvelopeBuild_c(SEXP bStarSEXP, SEXP ASEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP, SEXP nSEXP, SEXP sortgridSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(EnvelopeBuild_c(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, sortgrid));
    return rcpp_result_gen;
END_RCPP
}
// EnvelopeBuild_Ind_Normal_Gamma
List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar, NumericMatrix A, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, std::string family, std::string link, int Gridtype, int n, bool sortgrid);
RcppExport SEXP _glmbayes_EnvelopeBuild_Ind_Normal_Gamma(SEXP bStarSEXP, SEXP ASEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP, SEXP nSEXP, SEXP sortgridSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(EnvelopeBuild_Ind_Normal_Gamma(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, sortgrid));
    return rcpp_result_gen;
END_RCPP
}
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
// RSS
NumericVector RSS(NumericVector y, NumericMatrix x, NumericMatrix b, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_RSS(SEXP ySEXP, SEXP xSEXP, SEXP bSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(RSS(y, x, b, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// f2_gaussian
NumericVector f2_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_f2_gaussian(SEXP bSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(f2_gaussian(b, y, x, mu, P, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// Inv_f3_gaussian
arma::mat Inv_f3_gaussian(NumericMatrix cbars, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt);
RcppExport SEXP _glmbayes_Inv_f3_gaussian(SEXP cbarsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type cbars(cbarsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(Inv_f3_gaussian(cbars, y, x, mu, P, alpha, wt));
    return rcpp_result_gen;
END_RCPP
}
// rindep_norm_gamma_reg_std_cpp
Rcpp::List rindep_norm_gamma_reg_std_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rindep_norm_gamma_reg_std_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(rindep_norm_gamma_reg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// rindep_norm_gamma_reg_std_v2_cpp
Rcpp::List rindep_norm_gamma_reg_std_v2_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::List gamma_list, Rcpp::List UB_list, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rindep_norm_gamma_reg_std_v2_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP gamma_listSEXP, SEXP UB_listSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::List >::type gamma_list(gamma_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type UB_list(UB_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(rindep_norm_gamma_reg_std_v2_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// rindep_norm_gamma_reg_std_v3_cpp
Rcpp::List rindep_norm_gamma_reg_std_v3_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::List gamma_list, Rcpp::List UB_list, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rindep_norm_gamma_reg_std_v3_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP gamma_listSEXP, SEXP UB_listSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::List >::type gamma_list(gamma_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type UB_list(UB_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(rindep_norm_gamma_reg_std_v3_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// rindep_norm_gamma_reg_std_v4_cpp
Rcpp::List rindep_norm_gamma_reg_std_v4_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::List gamma_list, Rcpp::List UB_list, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rindep_norm_gamma_reg_std_v4_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP gamma_listSEXP, SEXP UB_listSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::List >::type gamma_list(gamma_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type UB_list(UB_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(rindep_norm_gamma_reg_std_v4_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// rindep_norm_gamma_reg_std_v5_cpp
Rcpp::List rindep_norm_gamma_reg_std_v5_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::List gamma_list, Rcpp::List UB_list, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rindep_norm_gamma_reg_std_v5_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP gamma_listSEXP, SEXP UB_listSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::List >::type gamma_list(gamma_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type UB_list(UB_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type progbar(progbarSEXP);
    rcpp_result_gen = Rcpp::wrap(rindep_norm_gamma_reg_std_v5_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// glmb_Standardize_Model
Rcpp::List glmb_Standardize_Model(NumericVector y, NumericMatrix x, NumericMatrix P, NumericMatrix bstar, NumericMatrix A1);
RcppExport SEXP _glmbayes_glmb_Standardize_Model(SEXP ySEXP, SEXP xSEXP, SEXP PSEXP, SEXP bstarSEXP, SEXP A1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bstar(bstarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A1(A1SEXP);
    rcpp_result_gen = Rcpp::wrap(glmb_Standardize_Model(y, x, P, bstar, A1));
    return rcpp_result_gen;
END_RCPP
}
// rnnorm_reg_std_cpp
Rcpp::List rnnorm_reg_std_cpp(int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, Function f2, Rcpp::List Envelope, Rcpp::CharacterVector family, Rcpp::CharacterVector link, int progbar);
RcppExport SEXP _glmbayes_rnnorm_reg_std_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP, SEXP f2SEXP, SEXP EnvelopeSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP progbarSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(rnnorm_reg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar));
    return rcpp_result_gen;
END_RCPP
}
// rnnorm_reg_cpp
Rcpp::List rnnorm_reg_cpp(int n, NumericVector y, NumericMatrix x, NumericVector mu, NumericMatrix P, NumericVector offset, NumericVector wt, double dispersion, Function f2, Function f3, NumericVector start, std::string family, std::string link, int Gridtype);
RcppExport SEXP _glmbayes_rnnorm_reg_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP offsetSEXP, SEXP wtSEXP, SEXP dispersionSEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP startSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< Function >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< Function >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type Gridtype(GridtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(rnnorm_reg_cpp(n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link, Gridtype));
    return rcpp_result_gen;
END_RCPP
}
// rnorm_reg_cpp
Rcpp::List rnorm_reg_cpp(int n, NumericVector y, NumericMatrix x, NumericVector mu, NumericMatrix P, NumericVector offset, NumericVector wt, double dispersion, Function f2, Function f3, NumericVector start, std::string family, std::string link, int Gridtype);
RcppExport SEXP _glmbayes_rnorm_reg_cpp(SEXP nSEXP, SEXP ySEXP, SEXP xSEXP, SEXP muSEXP, SEXP PSEXP, SEXP offsetSEXP, SEXP wtSEXP, SEXP dispersionSEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP startSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP GridtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< Function >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< Function >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type Gridtype(GridtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm_reg_cpp(n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link, Gridtype));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_glmbayes_EnvelopeBuild_c", (DL_FUNC) &_glmbayes_EnvelopeBuild_c, 13},
    {"_glmbayes_EnvelopeBuild_Ind_Normal_Gamma", (DL_FUNC) &_glmbayes_EnvelopeBuild_Ind_Normal_Gamma, 13},
    {"_glmbayes_Set_Grid", (DL_FUNC) &_glmbayes_Set_Grid, 3},
    {"_glmbayes_setlogP", (DL_FUNC) &_glmbayes_setlogP, 4},
    {"_glmbayes_RSS", (DL_FUNC) &_glmbayes_RSS, 5},
    {"_glmbayes_f2_gaussian", (DL_FUNC) &_glmbayes_f2_gaussian, 7},
    {"_glmbayes_Inv_f3_gaussian", (DL_FUNC) &_glmbayes_Inv_f3_gaussian, 7},
    {"_glmbayes_rindep_norm_gamma_reg_std_cpp", (DL_FUNC) &_glmbayes_rindep_norm_gamma_reg_std_cpp, 12},
    {"_glmbayes_rindep_norm_gamma_reg_std_v2_cpp", (DL_FUNC) &_glmbayes_rindep_norm_gamma_reg_std_v2_cpp, 14},
    {"_glmbayes_rindep_norm_gamma_reg_std_v3_cpp", (DL_FUNC) &_glmbayes_rindep_norm_gamma_reg_std_v3_cpp, 14},
    {"_glmbayes_rindep_norm_gamma_reg_std_v4_cpp", (DL_FUNC) &_glmbayes_rindep_norm_gamma_reg_std_v4_cpp, 14},
    {"_glmbayes_rindep_norm_gamma_reg_std_v5_cpp", (DL_FUNC) &_glmbayes_rindep_norm_gamma_reg_std_v5_cpp, 14},
    {"_glmbayes_glmb_Standardize_Model", (DL_FUNC) &_glmbayes_glmb_Standardize_Model, 5},
    {"_glmbayes_rnnorm_reg_std_cpp", (DL_FUNC) &_glmbayes_rnnorm_reg_std_cpp, 12},
    {"_glmbayes_rnnorm_reg_cpp", (DL_FUNC) &_glmbayes_rnnorm_reg_cpp, 14},
    {"_glmbayes_rnorm_reg_cpp", (DL_FUNC) &_glmbayes_rnorm_reg_cpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_glmbayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
