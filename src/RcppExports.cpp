// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Arrh_func_cpp
NumericVector Arrh_func_cpp(NumericVector A, NumericVector E, double R, double temp_K);
RcppExport SEXP _ABM_Arrh_func_cpp(SEXP ASEXP, SEXP ESEXP, SEXP RSEXP, SEXP temp_KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type temp_K(temp_KSEXP);
    rcpp_result_gen = Rcpp::wrap(Arrh_func_cpp(A, E, R, temp_K));
    return rcpp_result_gen;
END_RCPP
}
// CTM_cpp
Rcpp::NumericVector CTM_cpp(Rcpp::NumericVector tt, Rcpp::NumericVector top, Rcpp::NumericVector tmin, Rcpp::NumericVector tmax, Rcpp::NumericVector yopt);
RcppExport SEXP _ABM_CTM_cpp(SEXP ttSEXP, SEXP topSEXP, SEXP tminSEXP, SEXP tmaxSEXP, SEXP yoptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tt(ttSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type top(topSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type yopt(yoptSEXP);
    rcpp_result_gen = Rcpp::wrap(CTM_cpp(tt, top, tmin, tmax, yopt));
    return rcpp_result_gen;
END_RCPP
}
// grepl_cpp
LogicalVector grepl_cpp(std::string pattern, CharacterVector x);
RcppExport SEXP _ABM_grepl_cpp(SEXP patternSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pattern(patternSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(grepl_cpp(pattern, x));
    return rcpp_result_gen;
END_RCPP
}
// inhib_cpp
List inhib_cpp(bool pH_inhib_overrule, double pH, double NH3_frac, double HAC_frac, double H2S_frac, double TAN, double VFA, double sulfide, double slurry_mass, NumericVector pH_LL, NumericVector pH_UL, NumericVector ki_NH3_min, NumericVector ki_NH3_max, NumericVector ki_NH4_min, NumericVector ki_NH4_max, NumericVector ki_HAC, NumericVector ki_H2S_slope, NumericVector ki_H2S_int, NumericVector ki_H2S_min, NumericVector IC50_low);
RcppExport SEXP _ABM_inhib_cpp(SEXP pH_inhib_overruleSEXP, SEXP pHSEXP, SEXP NH3_fracSEXP, SEXP HAC_fracSEXP, SEXP H2S_fracSEXP, SEXP TANSEXP, SEXP VFASEXP, SEXP sulfideSEXP, SEXP slurry_massSEXP, SEXP pH_LLSEXP, SEXP pH_ULSEXP, SEXP ki_NH3_minSEXP, SEXP ki_NH3_maxSEXP, SEXP ki_NH4_minSEXP, SEXP ki_NH4_maxSEXP, SEXP ki_HACSEXP, SEXP ki_H2S_slopeSEXP, SEXP ki_H2S_intSEXP, SEXP ki_H2S_minSEXP, SEXP IC50_lowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type pH_inhib_overrule(pH_inhib_overruleSEXP);
    Rcpp::traits::input_parameter< double >::type pH(pHSEXP);
    Rcpp::traits::input_parameter< double >::type NH3_frac(NH3_fracSEXP);
    Rcpp::traits::input_parameter< double >::type HAC_frac(HAC_fracSEXP);
    Rcpp::traits::input_parameter< double >::type H2S_frac(H2S_fracSEXP);
    Rcpp::traits::input_parameter< double >::type TAN(TANSEXP);
    Rcpp::traits::input_parameter< double >::type VFA(VFASEXP);
    Rcpp::traits::input_parameter< double >::type sulfide(sulfideSEXP);
    Rcpp::traits::input_parameter< double >::type slurry_mass(slurry_massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pH_LL(pH_LLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pH_UL(pH_ULSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_NH3_min(ki_NH3_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_NH3_max(ki_NH3_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_NH4_min(ki_NH4_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_NH4_max(ki_NH4_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_HAC(ki_HACSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_H2S_slope(ki_H2S_slopeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_H2S_int(ki_H2S_intSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ki_H2S_min(ki_H2S_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type IC50_low(IC50_lowSEXP);
    rcpp_result_gen = Rcpp::wrap(inhib_cpp(pH_inhib_overrule, pH, NH3_frac, HAC_frac, H2S_frac, TAN, VFA, sulfide, slurry_mass, pH_LL, pH_UL, ki_NH3_min, ki_NH3_max, ki_NH4_min, ki_NH4_max, ki_HAC, ki_H2S_slope, ki_H2S_int, ki_H2S_min, IC50_low));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ABM_Arrh_func_cpp", (DL_FUNC) &_ABM_Arrh_func_cpp, 4},
    {"_ABM_CTM_cpp", (DL_FUNC) &_ABM_CTM_cpp, 5},
    {"_ABM_grepl_cpp", (DL_FUNC) &_ABM_grepl_cpp, 2},
    {"_ABM_inhib_cpp", (DL_FUNC) &_ABM_inhib_cpp, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_ABM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
