// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Arrh_func_cpp
NumericVector Arrh_func_cpp(NumericVector A, NumericVector E, double R, double temp_K, double scale_alpha_opt, double alpha_opt_scale_type, double alpha_opt_scale_CP);
RcppExport SEXP _ABM_Arrh_func_cpp(SEXP ASEXP, SEXP ESEXP, SEXP RSEXP, SEXP temp_KSEXP, SEXP scale_alpha_optSEXP, SEXP alpha_opt_scale_typeSEXP, SEXP alpha_opt_scale_CPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type temp_K(temp_KSEXP);
    Rcpp::traits::input_parameter< double >::type scale_alpha_opt(scale_alpha_optSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_opt_scale_type(alpha_opt_scale_typeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_opt_scale_CP(alpha_opt_scale_CPSEXP);
    rcpp_result_gen = Rcpp::wrap(Arrh_func_cpp(A, E, R, temp_K, scale_alpha_opt, alpha_opt_scale_type, alpha_opt_scale_CP));
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
// call_int
NumericVector call_int(Rcpp::Function temp_C_fun, double x);
RcppExport SEXP _ABM_call_int(SEXP temp_C_funSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type temp_C_fun(temp_C_funSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(call_int(temp_C_fun, x));
    return rcpp_result_gen;
END_RCPP
}
// combined_cpp
List combined_cpp(double g_NH4, bool pH_inhib_overrule, double pH, double pH_floor, double TAN, double VFA, double sulfide, double slurry_mass, NumericVector pH_LL, NumericVector pH_UL, NumericVector ki_NH3_min, NumericVector ki_NH3_max, NumericVector ki_NH4_min, NumericVector ki_NH4_max, NumericVector ki_HAC, NumericVector ki_H2S_slope, NumericVector ki_H2S_int, NumericVector ki_H2S_min, NumericVector IC50_low, double temp_K, double temp_C, double temp_standard, double area, double floor_area, double Cfat, double CPs, double CPf, double RFd, double starch, double VSd, bool resp, List kl, NumericVector qhat, NumericVector i_meth, NumericVector i_sr, NumericVector xa, NumericVector ks_coefficient, double scale_ks, double ks_SO4, double sulfate, double urea, double alpha_urea, double km_urea);
RcppExport SEXP _ABM_combined_cpp(SEXP g_NH4SEXP, SEXP pH_inhib_overruleSEXP, SEXP pHSEXP, SEXP pH_floorSEXP, SEXP TANSEXP, SEXP VFASEXP, SEXP sulfideSEXP, SEXP slurry_massSEXP, SEXP pH_LLSEXP, SEXP pH_ULSEXP, SEXP ki_NH3_minSEXP, SEXP ki_NH3_maxSEXP, SEXP ki_NH4_minSEXP, SEXP ki_NH4_maxSEXP, SEXP ki_HACSEXP, SEXP ki_H2S_slopeSEXP, SEXP ki_H2S_intSEXP, SEXP ki_H2S_minSEXP, SEXP IC50_lowSEXP, SEXP temp_KSEXP, SEXP temp_CSEXP, SEXP temp_standardSEXP, SEXP areaSEXP, SEXP floor_areaSEXP, SEXP CfatSEXP, SEXP CPsSEXP, SEXP CPfSEXP, SEXP RFdSEXP, SEXP starchSEXP, SEXP VSdSEXP, SEXP respSEXP, SEXP klSEXP, SEXP qhatSEXP, SEXP i_methSEXP, SEXP i_srSEXP, SEXP xaSEXP, SEXP ks_coefficientSEXP, SEXP scale_ksSEXP, SEXP ks_SO4SEXP, SEXP sulfateSEXP, SEXP ureaSEXP, SEXP alpha_ureaSEXP, SEXP km_ureaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type g_NH4(g_NH4SEXP);
    Rcpp::traits::input_parameter< bool >::type pH_inhib_overrule(pH_inhib_overruleSEXP);
    Rcpp::traits::input_parameter< double >::type pH(pHSEXP);
    Rcpp::traits::input_parameter< double >::type pH_floor(pH_floorSEXP);
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
    Rcpp::traits::input_parameter< double >::type temp_K(temp_KSEXP);
    Rcpp::traits::input_parameter< double >::type temp_C(temp_CSEXP);
    Rcpp::traits::input_parameter< double >::type temp_standard(temp_standardSEXP);
    Rcpp::traits::input_parameter< double >::type area(areaSEXP);
    Rcpp::traits::input_parameter< double >::type floor_area(floor_areaSEXP);
    Rcpp::traits::input_parameter< double >::type Cfat(CfatSEXP);
    Rcpp::traits::input_parameter< double >::type CPs(CPsSEXP);
    Rcpp::traits::input_parameter< double >::type CPf(CPfSEXP);
    Rcpp::traits::input_parameter< double >::type RFd(RFdSEXP);
    Rcpp::traits::input_parameter< double >::type starch(starchSEXP);
    Rcpp::traits::input_parameter< double >::type VSd(VSdSEXP);
    Rcpp::traits::input_parameter< bool >::type resp(respSEXP);
    Rcpp::traits::input_parameter< List >::type kl(klSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type i_meth(i_methSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type i_sr(i_srSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ks_coefficient(ks_coefficientSEXP);
    Rcpp::traits::input_parameter< double >::type scale_ks(scale_ksSEXP);
    Rcpp::traits::input_parameter< double >::type ks_SO4(ks_SO4SEXP);
    Rcpp::traits::input_parameter< double >::type sulfate(sulfateSEXP);
    Rcpp::traits::input_parameter< double >::type urea(ureaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_urea(alpha_ureaSEXP);
    Rcpp::traits::input_parameter< double >::type km_urea(km_ureaSEXP);
    rcpp_result_gen = Rcpp::wrap(combined_cpp(g_NH4, pH_inhib_overrule, pH, pH_floor, TAN, VFA, sulfide, slurry_mass, pH_LL, pH_UL, ki_NH3_min, ki_NH3_max, ki_NH4_min, ki_NH4_max, ki_HAC, ki_H2S_slope, ki_H2S_int, ki_H2S_min, IC50_low, temp_K, temp_C, temp_standard, area, floor_area, Cfat, CPs, CPf, RFd, starch, VSd, resp, kl, qhat, i_meth, i_sr, xa, ks_coefficient, scale_ks, ks_SO4, sulfate, urea, alpha_urea, km_urea));
    return rcpp_result_gen;
END_RCPP
}
// derivatives_cpp
NumericVector derivatives_cpp(NumericVector yield, NumericVector rut, double slurry_prod_rate, double decay_rate, NumericVector scale, NumericVector xa_fresh, List conc_fresh, List ferm, NumericVector alpha, NumericVector COD_conv, double rain, double evap, double area, double respiration, double sub_resp, double rut_urea, double NH3_emis_rate_pit, double NH3_emis_rate_floor, double N2O_emis_rate, NumericVector rutsr, double H2S_emis_rate, double R, double temp_K, NumericVector i_meth, NumericVector i_sr, double CO2_ferm_meth_sr, NumericVector xa, double slurry_mass, double xa_bac, double xa_aer, double xa_dead, double RFd, double iNDF, double ash, double VSd, double starch, double CPs, double CPf, double Cfat, double VFA, double urea, double TAN, double sulfate, double sulfide, NumericVector qhat);
RcppExport SEXP _ABM_derivatives_cpp(SEXP yieldSEXP, SEXP rutSEXP, SEXP slurry_prod_rateSEXP, SEXP decay_rateSEXP, SEXP scaleSEXP, SEXP xa_freshSEXP, SEXP conc_freshSEXP, SEXP fermSEXP, SEXP alphaSEXP, SEXP COD_convSEXP, SEXP rainSEXP, SEXP evapSEXP, SEXP areaSEXP, SEXP respirationSEXP, SEXP sub_respSEXP, SEXP rut_ureaSEXP, SEXP NH3_emis_rate_pitSEXP, SEXP NH3_emis_rate_floorSEXP, SEXP N2O_emis_rateSEXP, SEXP rutsrSEXP, SEXP H2S_emis_rateSEXP, SEXP RSEXP, SEXP temp_KSEXP, SEXP i_methSEXP, SEXP i_srSEXP, SEXP CO2_ferm_meth_srSEXP, SEXP xaSEXP, SEXP slurry_massSEXP, SEXP xa_bacSEXP, SEXP xa_aerSEXP, SEXP xa_deadSEXP, SEXP RFdSEXP, SEXP iNDFSEXP, SEXP ashSEXP, SEXP VSdSEXP, SEXP starchSEXP, SEXP CPsSEXP, SEXP CPfSEXP, SEXP CfatSEXP, SEXP VFASEXP, SEXP ureaSEXP, SEXP TANSEXP, SEXP sulfateSEXP, SEXP sulfideSEXP, SEXP qhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yield(yieldSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rut(rutSEXP);
    Rcpp::traits::input_parameter< double >::type slurry_prod_rate(slurry_prod_rateSEXP);
    Rcpp::traits::input_parameter< double >::type decay_rate(decay_rateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xa_fresh(xa_freshSEXP);
    Rcpp::traits::input_parameter< List >::type conc_fresh(conc_freshSEXP);
    Rcpp::traits::input_parameter< List >::type ferm(fermSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type COD_conv(COD_convSEXP);
    Rcpp::traits::input_parameter< double >::type rain(rainSEXP);
    Rcpp::traits::input_parameter< double >::type evap(evapSEXP);
    Rcpp::traits::input_parameter< double >::type area(areaSEXP);
    Rcpp::traits::input_parameter< double >::type respiration(respirationSEXP);
    Rcpp::traits::input_parameter< double >::type sub_resp(sub_respSEXP);
    Rcpp::traits::input_parameter< double >::type rut_urea(rut_ureaSEXP);
    Rcpp::traits::input_parameter< double >::type NH3_emis_rate_pit(NH3_emis_rate_pitSEXP);
    Rcpp::traits::input_parameter< double >::type NH3_emis_rate_floor(NH3_emis_rate_floorSEXP);
    Rcpp::traits::input_parameter< double >::type N2O_emis_rate(N2O_emis_rateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rutsr(rutsrSEXP);
    Rcpp::traits::input_parameter< double >::type H2S_emis_rate(H2S_emis_rateSEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type temp_K(temp_KSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type i_meth(i_methSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type i_sr(i_srSEXP);
    Rcpp::traits::input_parameter< double >::type CO2_ferm_meth_sr(CO2_ferm_meth_srSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< double >::type slurry_mass(slurry_massSEXP);
    Rcpp::traits::input_parameter< double >::type xa_bac(xa_bacSEXP);
    Rcpp::traits::input_parameter< double >::type xa_aer(xa_aerSEXP);
    Rcpp::traits::input_parameter< double >::type xa_dead(xa_deadSEXP);
    Rcpp::traits::input_parameter< double >::type RFd(RFdSEXP);
    Rcpp::traits::input_parameter< double >::type iNDF(iNDFSEXP);
    Rcpp::traits::input_parameter< double >::type ash(ashSEXP);
    Rcpp::traits::input_parameter< double >::type VSd(VSdSEXP);
    Rcpp::traits::input_parameter< double >::type starch(starchSEXP);
    Rcpp::traits::input_parameter< double >::type CPs(CPsSEXP);
    Rcpp::traits::input_parameter< double >::type CPf(CPfSEXP);
    Rcpp::traits::input_parameter< double >::type Cfat(CfatSEXP);
    Rcpp::traits::input_parameter< double >::type VFA(VFASEXP);
    Rcpp::traits::input_parameter< double >::type urea(ureaSEXP);
    Rcpp::traits::input_parameter< double >::type TAN(TANSEXP);
    Rcpp::traits::input_parameter< double >::type sulfate(sulfateSEXP);
    Rcpp::traits::input_parameter< double >::type sulfide(sulfideSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    rcpp_result_gen = Rcpp::wrap(derivatives_cpp(yield, rut, slurry_prod_rate, decay_rate, scale, xa_fresh, conc_fresh, ferm, alpha, COD_conv, rain, evap, area, respiration, sub_resp, rut_urea, NH3_emis_rate_pit, NH3_emis_rate_floor, N2O_emis_rate, rutsr, H2S_emis_rate, R, temp_K, i_meth, i_sr, CO2_ferm_meth_sr, xa, slurry_mass, xa_bac, xa_aer, xa_dead, RFd, iNDF, ash, VSd, starch, CPs, CPf, Cfat, VFA, urea, TAN, sulfate, sulfide, qhat));
    return rcpp_result_gen;
END_RCPP
}
// extract_xa_cpp
NumericVector extract_xa_cpp(NumericVector y, int n_mic);
RcppExport SEXP _ABM_extract_xa_cpp(SEXP ySEXP, SEXP n_micSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n_mic(n_micSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_xa_cpp(y, n_mic));
    return rcpp_result_gen;
END_RCPP
}
// rates_cpp
List rates_cpp(double t, NumericVector y, List parms, Rcpp::Function temp_C_fun, Rcpp::Function pH_fun, Rcpp::Function SO4_inhibition_fun, List conc_fresh_fun, NumericVector xa_fresh_fun, Rcpp::Function CTM_cpp);
RcppExport SEXP _ABM_rates_cpp(SEXP tSEXP, SEXP ySEXP, SEXP parmsSEXP, SEXP temp_C_funSEXP, SEXP pH_funSEXP, SEXP SO4_inhibition_funSEXP, SEXP conc_fresh_funSEXP, SEXP xa_fresh_funSEXP, SEXP CTM_cppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type temp_C_fun(temp_C_funSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pH_fun(pH_funSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type SO4_inhibition_fun(SO4_inhibition_funSEXP);
    Rcpp::traits::input_parameter< List >::type conc_fresh_fun(conc_fresh_funSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xa_fresh_fun(xa_fresh_funSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type CTM_cpp(CTM_cppSEXP);
    rcpp_result_gen = Rcpp::wrap(rates_cpp(t, y, parms, temp_C_fun, pH_fun, SO4_inhibition_fun, conc_fresh_fun, xa_fresh_fun, CTM_cpp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ABM_Arrh_func_cpp", (DL_FUNC) &_ABM_Arrh_func_cpp, 7},
    {"_ABM_CTM_cpp", (DL_FUNC) &_ABM_CTM_cpp, 5},
    {"_ABM_call_int", (DL_FUNC) &_ABM_call_int, 2},
    {"_ABM_combined_cpp", (DL_FUNC) &_ABM_combined_cpp, 43},
    {"_ABM_derivatives_cpp", (DL_FUNC) &_ABM_derivatives_cpp, 45},
    {"_ABM_extract_xa_cpp", (DL_FUNC) &_ABM_extract_xa_cpp, 2},
    {"_ABM_rates_cpp", (DL_FUNC) &_ABM_rates_cpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_ABM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
