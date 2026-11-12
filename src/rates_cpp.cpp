#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include "ctm.h"

using namespace Rcpp;

// Derivatives for the ABM model, called by deSolve::lsoda() through abm_regular() and abm_variable()
//
// Notes on speed: this function is called very many times (~100 000 times for a 10 year simulation)
// so it avoids 1) calls back to R, 2) allocation of R vectors (except for the returned ones),
// and 3) name lookups in the parameter list (parameters are accessed by position through p_idx,
// see pars_indices.R). Time-variable inputs are interpolated here from tables made by makeInterpTab().

// Interpolation equivalent to stats::approxfun(x, y, method, rule = 2, f = 0)
// method 0 = linear, 1 = constant (f = 0, i.e., value of left knot)
// x must be sorted and without duplicates (ensured by makeInterpTab())
static inline double interp1(double xout, const double *x, const double *y, int n, int method) {
  if (n == 1 || xout <= x[0]) return y[0];
  if (xout >= x[n - 1]) return y[n - 1];
  // First knot > xout
  int j = std::upper_bound(x, x + n, xout) - x;
  int i = j - 1;
  if (method == 1) return y[i];
  return y[i] + (y[j] - y[i]) * (xout - x[i]) / (x[j] - x[i]);
}

// Interpolate column `col` of a table created by makeInterpTab() in R:
// list(time = numeric, y = numeric matrix (n time x k columns), method = integer)
static inline double interpTab(SEXP tab, double xout, int col = 0) {
  SEXP x = VECTOR_ELT(tab, 0);
  SEXP y = VECTOR_ELT(tab, 1);
  int method = Rf_asInteger(VECTOR_ELT(tab, 2));
  int n = Rf_length(x);
  // Matrix is column-major, so column col starts at col * n
  return interp1(xout, REAL(x), REAL(y) + (R_xlen_t)col * n, n, method);
}

static inline int interpTabNcol(SEXP tab) {
  return Rf_ncols(VECTOR_ELT(tab, 1));
}

// pH from sulfate concentration for acidified pig slurry (same as H2SO4_titrat(), class_anim = "pig")
static double pH_from_SO4_pig(double conc_SO4) {
  static const double cSO4[14] = {0, 0.237188238, 0.474376476, 0.711564714, 0.948752952,
                                  1.18594119, 1.423129427, 1.660317665, 1.897505903,
                                  2.134694141, 2.371882379, 2.609070617, 2.846258855,
                                  3.083447093};
  static const double pHv[14] = {7, 6.6375, 6.435, 6.2925, 6.1475, 5.9975, 5.845,
                                 5.67, 5.495, 5.325, 5.1575, 5, 4.875, 4.755};
  return interp1(conc_SO4, cSO4, pHv, 14, 0);
}

// rng = false: no need to save/restore R's random number state (this is slow for many calls)
// [[Rcpp::export(rng = false)]]
List rates_cpp(double t, NumericVector y_in, List parms){

  // Parameter access by position (see pars_indices.R)
  NumericVector p_idx = parms["p_idx"];
  auto P = [&](int i) -> SEXP { return VECTOR_ELT(parms, (R_xlen_t)p_idx[i]); };
  auto D = [&](int i) -> double { return Rf_asReal(P(i)); };
  // Pointer to numeric vector (all vector parameters are REALSXP, checked in R)
  auto V = [&](int i) -> const double* { return REAL(P(i)); };

  int n_mic = Rf_asInteger(P(2));

  // Local copy of state variables, with negative or zero values set to a small positive value
  const int ny = y_in.size();
  std::vector<double> y(y_in.begin(), y_in.end());
  for (int i = 0; i < ny; ++i) {
    if (y[i] < 1E-10) y[i] = 1E-10;
  }

  double slurry_prod_rate = D(0);
  double t_run = D(1);

  // Extract state variable values from y argument
  const double *xa = &y[0];
  double slurry_mass = y[n_mic];
  double xa_dead = y[n_mic+1];
  double RFd = y[n_mic+2];
  double VSd = y[n_mic+5];
  double starch = y[n_mic+6];
  double CPs = y[n_mic+7];
  double CPf = y[n_mic+8];
  double Cfat = y[n_mic+9];
  double VFA = y[n_mic+10];
  double urea = y[n_mic+11];
  double TAN = y[n_mic+12];
  double sulfate = y[n_mic+13];
  double sulfide = y[n_mic+14];

  double time_sum = t + t_run; // Time since start of complete simulation

  // Adjust slurry_prod_rate if grazing
  double t_max_graze = D(58);
  if (t_max_graze != 0.0){
    double t_year = time_sum - floor(time_sum/365) * 365;
    double t_min_graze = D(57);
    double slurry_graze_reduc = D(59);
    if(t_year >= t_min_graze && t_year <= t_max_graze){
      slurry_prod_rate *= slurry_graze_reduc;
    }
  }

  // pH and temperature from interpolation tables (constant input = table with 1 row)
  double pH = interpTab(P(61), time_sum);
  // pH = -1 is code for calculation from sulfate concentration (pH = 'calc')
  if (pH == -1.0) {
    pH = pH_from_SO4_pig(sulfate / slurry_mass);
  }

  double temp_C = interpTab(P(60), time_sum);
  if (temp_C <= 0.0) temp_C = 0.00001;
  double temp_K = temp_C + 273.15;
  double temp_standard = D(3);

  // Arrhenius parms for hydrolysis
  const double *A = V(4);
  const double *E = V(5);
  int n_alpha = Rf_length(P(4));
  const double *scale = V(6);
  double scale_ks_coefficient = scale[0];
  double scale_qhat_opt = scale[1];
  double scale_xa_fresh = scale[2];
  double scale_yield = scale[3];
  double scale_alpha_opt = scale[4];

  double R = D(7);
  double alpha_opt_scale_type = D(8);
  double alpha_opt_scale_CP = D(9);

  // Hydrolysis rates with Arrhenius function
  // Order: xa_dead, starch, Cfat, CPs, CPf, RFd, VSd, urea
  std::vector<double> alpha(n_alpha);
  for (int i = 0; i < n_alpha; ++i) {
    alpha[i] = A[i] * exp(-E[i] / (R * temp_K));
  }
  for (int i = 0; i <= 6; ++i) alpha[i] = alpha[i] * scale_alpha_opt * alpha_opt_scale_type;
  for (int i = 3; i <= 4; ++i) alpha[i] = alpha[i] * alpha_opt_scale_CP;

  double alpha_xa_dead = alpha[0];
  double alpha_starch = alpha[1];
  double alpha_Cfat = alpha[2];
  double alpha_CPs = alpha[3];
  double alpha_CPf = alpha[4];
  double alpha_RFd = alpha[5];
  double alpha_VSd = alpha[6];
  double alpha_urea = alpha[7];

  const double *T_opt = V(10);
  const double *T_min = V(11);
  const double *T_max = V(12);
  const double *qhat_opt = V(13);

  // Microbial substrate utilization rate
  std::vector<double> qhat(n_mic);
  for (int i = 0; i < n_mic; ++i) {
    qhat[i] = scale_qhat_opt * ctm1(temp_K, T_opt[i], T_min[i], T_max[i], qhat_opt[i]);
  }

  // Decay of microorganisms: group-specific rate (decay_rate) at 313 K, faster decay at higher temp up until 313 K, 
  // above which decay is constant at decay_rate
  const double *decay_rate_opt = V(14);
  std::vector<double> decay_rate(n_mic);
  for (int i = 0; i < n_mic; ++i) {
    decay_rate[i] = (temp_K > 313) ? decay_rate_opt[i] : ctm1(temp_K, 313, 273, 325, decay_rate_opt[i]);
  }

  double g_NH4 = D(15);
  bool pH_inhib_overrule = Rf_asLogical(P(16));
  double pH_floor = D(17);
  const double *pH_LL = V(18);
  const double *pH_UL = V(19);
  const double *ki_NH3_min = V(20);
  const double *ki_NH3_max = V(21);
  const double *ki_NH4_min = V(22);
  const double *ki_NH4_max = V(23);
  const double *ki_HAC = V(24);
  const double *ki_H2S_slope = V(25);
  const double *ki_H2S_int = V(26);
  const double *ki_H2S_min = V(27);
  // V(28) = IC50_low: no longer used (IC50 is continuous at pH 6.8), slot kept so indices stay fixed
  double area = D(29);
  double floor_area = D(30);
  bool resp = Rf_asLogical(P(31));
  const double *kl = V(32);
  const double *i_meth = V(33);
  int n_meth = Rf_length(P(33));
  const double *i_sr = V(34);
  int n_sr = Rf_length(P(34));
  const double *ks_coefficient = V(35);
  double ks_SO4 = D(36);
  double km_urea = D(37);

  double log_HAC = -4.8288 + 21.42/temp_K;
  double log_NH3 = - 0.09046 - 2729.31/temp_K;
  double log_H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K);

  // Compute the fractions
  double HAC_frac = 1 - (1 / (1 + std::pow(10, -log_HAC - pH)));
  double NH3_frac = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH));
  double NH3_frac_floor = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH_floor));
  double H2S_frac = 1 - (1 / (1 + std::pow(10, -log_H2S - pH)));

  // Inhibition (all 0 if pH_inhib_overrule, except cum_inhib)
  std::vector<double> NH3_inhib(n_mic), NH4_inhib(n_mic), HAC_inhib(n_mic), H2S_inhib(n_mic), cum_inhib(n_mic);

  if (pH_inhib_overrule) {
    // Compute pH inhibition and set cumulative inhibition = pH_inhib
    for (int i = 0; i < n_mic; ++i) {
      cum_inhib[i] = (1 + 2 * pow(10, 0.5 * (pH_LL[i] - pH_UL[i]))) /
        (1 + pow(10, pH - pH_UL[i]) + pow(10, pH_LL[i] - pH));
    }
  } else {
    double NH3_conc = NH3_frac * TAN / slurry_mass;
    double NH4_conc = (1 - NH3_frac) * TAN / slurry_mass;
    double HAC_conc = HAC_frac * VFA / slurry_mass;
    double x = H2S_frac * sulfide / slurry_mass;

    for (int i = 0; i < n_mic; ++i) {
      // NH3 inhibition
      NH3_inhib[i] = (NH3_conc <= ki_NH3_min[i]) ? 1.0 :
      exp(-2.77259 * pow((NH3_conc - ki_NH3_min[i]) / (ki_NH3_max[i] - ki_NH3_min[i]), 2));

      // NH4 inhibition
      NH4_inhib[i] = (NH4_conc <= ki_NH4_min[i]) ? 1.0 :
        exp(-2.77259 * pow((NH4_conc - ki_NH4_min[i]) / (ki_NH4_max[i] - ki_NH4_min[i]), 2));

      // HAC inhibition
      HAC_inhib[i] = (HAC_conc >= 0.05) ?
      ((2 - ki_HAC[i] / (ki_HAC[i] + 0.05)) * ki_HAC[i] / (ki_HAC[i] + HAC_conc)) :
        1.0;

      // H2S inhibition
      // Linear in free H2S (x, g/kg): 1 at x = H2S_frac * ki_H2S_min, 0.5 at x = IC50
      // IC50 is held at its pH 6.8 value below pH 6.8 (continuous in pH)
      double IC50 = ki_H2S_slope[i] * std::max(pH, 6.8) + ki_H2S_int[i];
      double a = -0.5 / (IC50 - (H2S_frac * ki_H2S_min[i]));
      double b = 1 - a * H2S_frac * ki_H2S_min[i];

      H2S_inhib[i] = a * x + b;
      H2S_inhib[i] = std::max(0.0, std::min(1.0, H2S_inhib[i]));  // Clamp between 0 and 1

      // Compute cumulative inhibition
      cum_inhib[i] = HAC_inhib[i] * NH3_inhib[i] * NH4_inhib[i] * H2S_inhib[i];
    }
  }

  // Henry's constant for NH3
  double H_NH3 = 1431 * pow(1.053, (293 - temp_K));

  // Henry's constant for oxygen
  double kH_oxygen = 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000;

  // Mass transfer coefficients (order NH3, NH3_floor, H2S)
  double kl_NH3 = kl[0];
  double kl_NH3_floor = kl[1];
  double kl_H2S = kl[2];

  // NH3 emission rates
  double NH3_emis_rate_floor = kl_NH3_floor * floor_area * ((NH3_frac_floor * TAN) / slurry_mass) * 1000 / H_NH3;
  double NH3_emis_rate_pit = kl_NH3 * area * ((NH3_frac * TAN) / slurry_mass) * 1000 / H_NH3;

  // H2S emission rate
  double H2S_emis_rate = kl_H2S * area * ((H2S_frac * sulfide) / slurry_mass) * 1000;

  // Respiration
  double respiration = 0.0;
  double sub_resp = Cfat + CPs + CPf + RFd + starch + VSd;

  if (resp && (slurry_mass / area >= 1)) {
    double kl_oxygen = exp(0.6158816 + 0.09205127 * temp_C);
    respiration = kl_oxygen * area * (kH_oxygen * 0.208) * (sub_resp / slurry_mass) / 100;
  }

  // ks is temp dependent
  std::vector<double> ks(n_mic);
  for (int i = 0; i < n_mic; ++i) {
    ks[i] = ks_coefficient[i] * (0.8157 * exp(-0.063 * temp_C));
  }

  // VFA consumption rates (g/d)
  std::vector<double> rut(n_mic, 0.0);

  // By sulfate reducers, affected by inhibition terms
  for (int i = 0; i < n_sr; ++i) {
    int idx_sr = (int) i_sr[i];
    rut[idx_sr] = ((qhat[idx_sr] * VFA / slurry_mass * xa[idx_sr] / slurry_mass / (scale_ks_coefficient * ks[idx_sr] + VFA / slurry_mass)) * slurry_mass *
      (sulfate / slurry_mass) / (ks_SO4 + sulfate / slurry_mass)) * cum_inhib[idx_sr];
  }

  // By methanogens, affected by inhibition terms
  for(int i = 0; i < n_meth; ++i) {
    int idx = (int) i_meth[i];
    rut[idx] = ((qhat[idx] * VFA / (slurry_mass) * xa[idx] / (slurry_mass)) / (scale_ks_coefficient * ks[idx] + VFA / (slurry_mass)) *
      (slurry_mass)) * cum_inhib[idx];
  }

  // Stop if any rut are below 0
  for (int i = 0; i < n_mic; ++i) {
    if (rut[i] < 0) Rcpp::stop("In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)");
  }

  // Urea hydrolysis by Michaelis Menten
  double rut_urea = (alpha_urea * urea) / (km_urea + urea/slurry_mass);

  // Fresh slurry composition, from interpolation table (constant input = table with 1 row)
  // Column order is fixed in R (see conc_fresh_names in internal_helpers.R)
  SEXP cf_tab = P(62);
  int n_cf = interpTabNcol(cf_tab);
  double conc_fresh[14];
  for (int i = 0; i < 14; ++i) conc_fresh[i] = interpTab(cf_tab, time_sum, i);

  double conc_fresh_sulfide = conc_fresh[0];
  double conc_fresh_urea = conc_fresh[1];
  double conc_fresh_sulfate = conc_fresh[2];
  double conc_fresh_TAN = conc_fresh[3];
  double conc_fresh_starch = conc_fresh[4];
  double conc_fresh_VFA = conc_fresh[5];
  double conc_fresh_xa_dead = conc_fresh[6];
  double conc_fresh_Cfat = conc_fresh[7];
  double conc_fresh_CPs = conc_fresh[8];
  double conc_fresh_CPf = conc_fresh[9];
  double conc_fresh_RFd = conc_fresh[10];
  double conc_fresh_iNDF = conc_fresh[11];
  double conc_fresh_VSd = conc_fresh[12];
  double conc_fresh_ash = conc_fresh[13];

  // Fermentation stoichiometry
  double mol_carb = (alpha_RFd * RFd + alpha_starch * starch) * 0.005208333;
  double mol_pro = (alpha_CPf * CPf + alpha_CPs * CPs) * 0.00748503;
  double mol_lip = alpha_Cfat * Cfat * 0.0004194631;
  double mol_OM = alpha_VSd * VSd * 0.001833952;

  const double *carb = V(39);
  int n_carb = Rf_length(P(39));
  const double *pro = V(40);
  const double *lip = V(41);
  const double *OM = V(55);
  const double *ace = V(42);
  const double *hyd = V(43);
  const double *ace_sr = V(44);
  const double *hyd_sr = V(45);

  std::vector<double> ferm(n_carb);
  for(int i = 0; i < n_carb; ++i) {
    ferm[i] = mol_carb * carb[i] + mol_pro * pro[i] + mol_lip * lip[i] + mol_OM * OM[i];
  }

  double C2H4O2 = ferm[6];
  double H2 = ferm[7];

  double VFA_H2_ferm = C2H4O2 * 64 + H2 * 16;
  double TAN_min_ferm = ferm[3] * 14.007;

  double COD_conv_meth_CO2 = -(ace[1] * 64 + hyd[0] * 16)/((ace[2] + hyd[2]) * 44.01);
  double COD_conv_sr_CO2 = -(ace_sr[1] * 64 + hyd_sr[0] * 16)/((ace_sr[3]) * 44.01);
  double CO2_ferm = ferm[8] * 44.01;

  // Respiration
  double mol_carb_resp = respiration * (RFd + starch)/sub_resp * 0.005208333;
  double mol_pro_resp = respiration * (CPs + CPf)/sub_resp * 0.00748503;
  double mol_lip_resp = respiration * Cfat/sub_resp * 0.0004194631;
  double mol_OM_resp = respiration * VSd/sub_resp * 0.001833952;

  const double *carb_resp = V(46);
  int n_carb_resp = Rf_length(P(46));
  const double *pro_resp = V(47);
  const double *lip_resp = V(48);
  const double *OM_resp = V(56);

  std::vector<double> respir(n_carb_resp);
  for(int i = 0; i < n_carb_resp; ++i) {
    respir[i] = mol_carb_resp * carb_resp[i] + mol_pro_resp * pro_resp[i] + mol_lip_resp * lip_resp[i] + mol_OM_resp * OM_resp[i];
  }

  double TAN_min_resp = respir[3] * 14.007;
  double CO2_resp = (respir[8] + respir[4]) * 44.01;

  // COD and C conversion factors (order: see COD_conv_names in internal_helpers.R)
  const double *COD_conv = V(49);
  double COD_conv_S = COD_conv[8];
  double COD_conv_CH4 = COD_conv[0];
  double COD_conv_CO2_ureo = COD_conv[12];
  double COD_conv_CP_N = COD_conv[13];
  double COD_conv_C_RFd = COD_conv[15];
  double COD_conv_C_iNDF = COD_conv[16];
  double COD_conv_C_starch = COD_conv[18];
  double COD_conv_C_VFA = COD_conv[21];
  double COD_conv_C_xa = COD_conv[14];
  double COD_conv_C_Cfat = COD_conv[19];
  double COD_conv_C_CP = COD_conv[20];
  double COD_conv_C_VSd = COD_conv[22];
  double COD_conv_C_N_urea = COD_conv[23];

  double sum_rut = std::accumulate(rut.begin(), rut.end(), 0.0);
  double sum_rutmeth = 0.0;
  for(int i = 0; i < n_meth; ++i) {
    sum_rutmeth += rut[(int) i_meth[i]];
  }
  double sum_rutsr = 0.0;
  for(int i = 0; i < n_sr; ++i) {
    sum_rutsr += rut[(int) i_sr[i]];
  }

  // CO2 produced from combined fermentation, methanogenesis and sulfate reduction
  double CO2_ferm_meth_sr = CO2_ferm + sum_rutmeth/COD_conv_meth_CO2 + sum_rutsr/COD_conv_sr_CO2;

  // Microbial biomass in fresh slurry, from interpolation table (columns in order of microbial groups)
  SEXP xaf_tab = P(63);
  std::vector<double> xa_fresh(n_mic);
  for (int i = 0; i < n_mic; ++i) xa_fresh[i] = interpTab(xaf_tab, time_sum, i);
  double sum_xa_fresh = std::accumulate(xa_fresh.begin(), xa_fresh.end(), 0.0);
  double sum_decay = 0.0;
  for (int i = 0; i < n_mic; ++i) sum_decay += decay_rate[i] * xa[i];

  const double *yield = V(51);
  double rain = D(52);
  double evap = D(53);
  double N2O_emis_rate = D(54);

  // Mass balance diagnostics: element consumed - element produced by each process (g/d)
  // Positive = element destroyed by process stoichiometry, negative = created
  // Computed from process rates (not from derivatives) so the balance also checks the derivatives
  const double C_CH4 = 12.011 / 16.043;
  const double C_CO2 = 12.011 / 44.01;
  double h_RFd = alpha_RFd * RFd;
  double h_starch = alpha_starch * starch;
  double h_CP = alpha_CPs * CPs + alpha_CPf * CPf;
  double h_Cfat = alpha_Cfat * Cfat;
  double h_VSd = alpha_VSd * VSd;
  double h_xa_dead = alpha_xa_dead * xa_dead;
  double sum_growth = 0.0;
  double sum_growth_C = 0.0;
  for (int i = 0; i < n_mic; ++i) sum_growth += scale_yield * yield[i] * rut[i];
  sum_growth_C = sum_growth / COD_conv_C_xa;

  // Hydrolysis + fermentation: substrate -> VFA (incl. H2) + CO2 + TAN
  double COD_gap_ferm = h_RFd + h_starch + h_CP + h_Cfat + h_VSd - VFA_H2_ferm;
  double C_gap_ferm = h_RFd / COD_conv_C_RFd + h_starch / COD_conv_C_starch + h_CP / COD_conv_C_CP +
    h_Cfat / COD_conv_C_Cfat + h_VSd / COD_conv_C_VSd - VFA_H2_ferm / COD_conv_C_VFA - CO2_ferm * C_CO2;
  double N_gap_ferm = h_CP / COD_conv_CP_N - TAN_min_ferm;
  // Microbial uptake of VFA: VFA -> CH4 + CO2 + sulfide (COD) + biomass
  double COD_gap_growth = -sum_growth;
  double C_gap_uptake = sum_rut / COD_conv_C_VFA - sum_rutmeth / COD_conv_CH4 * C_CH4 -
    (sum_rutmeth / COD_conv_meth_CO2 + sum_rutsr / COD_conv_sr_CO2) * C_CO2 - sum_growth_C;
  // Hydrolysis of dead biomass: xa_dead -> VFA
  double C_gap_xa_dead = h_xa_dead / COD_conv_C_xa - h_xa_dead / COD_conv_C_VFA;
  // Respiration: substrate -> CO2 + TAN (COD is tracked as COD_conv_cum_respir)
  double C_gap_respir = 0.0;
  double N_gap_respir = 0.0;
  if (respiration > 0.0) {
    C_gap_respir = respiration / sub_resp * (RFd / COD_conv_C_RFd + starch / COD_conv_C_starch + (CPs + CPf) / COD_conv_C_CP +
      Cfat / COD_conv_C_Cfat + VSd / COD_conv_C_VSd) - CO2_resp * C_CO2;
    N_gap_respir = respiration / sub_resp * (CPs + CPf) / COD_conv_CP_N - TAN_min_resp;
  }
  // Urea hydrolysis: urea -> TAN + CO2
  double C_gap_urea = rut_urea / COD_conv_C_N_urea - rut_urea / COD_conv_CO2_ureo * C_CO2;

  // Derivatives
  int nn = n_mic + 39;
  NumericVector derivatives(nn);

  for (int i = 0; i < n_mic; i++) {
    derivatives[i] = scale_yield * yield[i] * rut[i] + scale_xa_fresh * xa_fresh[i] * slurry_prod_rate - decay_rate[i] * xa[i];
  }

  derivatives[n_mic] = slurry_prod_rate + (rain - evap) * area;
  derivatives[n_mic+1] = slurry_prod_rate * conc_fresh_xa_dead - alpha_xa_dead * xa_dead + sum_decay;
  derivatives[n_mic+2] = slurry_prod_rate * conc_fresh_RFd - alpha_RFd * RFd - respiration * RFd / sub_resp;
  derivatives[n_mic+3] = slurry_prod_rate * conc_fresh_iNDF;
  derivatives[n_mic+4] = slurry_prod_rate * conc_fresh_ash;
  derivatives[n_mic+5] = slurry_prod_rate * conc_fresh_VSd - alpha_VSd * VSd - respiration * VSd / sub_resp;
  derivatives[n_mic+6] = slurry_prod_rate * conc_fresh_starch - alpha_starch * starch - respiration * starch / sub_resp;
  derivatives[n_mic+7] = slurry_prod_rate * conc_fresh_CPs - alpha_CPs * CPs - respiration * CPs / sub_resp;
  derivatives[n_mic+8] = slurry_prod_rate * conc_fresh_CPf - alpha_CPf * CPf - respiration * CPf / sub_resp;
  derivatives[n_mic+9] = slurry_prod_rate * conc_fresh_Cfat - alpha_Cfat * Cfat - respiration * Cfat / sub_resp;
  derivatives[n_mic+10] = alpha_xa_dead * xa_dead + VFA_H2_ferm - sum_rut + slurry_prod_rate * conc_fresh_VFA;
  derivatives[n_mic+11] = slurry_prod_rate * conc_fresh_urea - rut_urea;
  derivatives[n_mic+12] = slurry_prod_rate * conc_fresh_TAN + rut_urea + TAN_min_ferm + TAN_min_resp - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate;
  derivatives[n_mic+13] = slurry_prod_rate * conc_fresh_sulfate - sum_rutsr / COD_conv_S;
  derivatives[n_mic+14] = slurry_prod_rate * conc_fresh_sulfide + sum_rutsr / COD_conv_S - H2S_emis_rate;
  derivatives[n_mic+15] = NH3_emis_rate_pit + NH3_emis_rate_floor;
  derivatives[n_mic+16] = N2O_emis_rate;
  derivatives[n_mic+17] = sum_rutmeth / COD_conv_CH4;
  derivatives[n_mic+18] = CO2_ferm_meth_sr + CO2_resp + rut_urea / COD_conv_CO2_ureo;
  derivatives[n_mic+19] = sum_rut + respiration;
  derivatives[n_mic+20] = sum_rutmeth;
  derivatives[n_mic+21] = respiration;
  derivatives[n_mic+22] = sum_rutsr;
  derivatives[n_mic+23] = slurry_prod_rate * (conc_fresh_starch + conc_fresh_VFA + conc_fresh_xa_dead + conc_fresh_Cfat + conc_fresh_CPs + conc_fresh_CPf + conc_fresh_RFd + conc_fresh_iNDF + conc_fresh_VSd + (sum_xa_fresh * scale_xa_fresh));
  derivatives[n_mic+24] = slurry_prod_rate * (conc_fresh_starch/COD_conv_C_starch + conc_fresh_RFd/COD_conv_C_RFd + conc_fresh_iNDF/COD_conv_C_iNDF + conc_fresh_VFA/COD_conv_C_VFA + conc_fresh_xa_dead/COD_conv_C_xa +
    conc_fresh_Cfat/COD_conv_C_Cfat + (conc_fresh_CPs + conc_fresh_CPf)/COD_conv_C_CP +
    conc_fresh_VSd/COD_conv_C_VSd + conc_fresh_urea/COD_conv_C_N_urea + (sum_xa_fresh * scale_xa_fresh)/(COD_conv_C_xa));
  derivatives[n_mic+25] =  slurry_prod_rate * ((conc_fresh_CPs + conc_fresh_CPf)/COD_conv_CP_N + conc_fresh_TAN + conc_fresh_urea);
  derivatives[n_mic+26] = slurry_prod_rate;
  // Mass balance: n_mic+27 to n_mic+29 are COD_eff_cum, C_eff_cum, N_eff_cum (changed only at emptying, derivative 0)
  derivatives[n_mic+30] = COD_gap_ferm;
  derivatives[n_mic+31] = COD_gap_growth;
  derivatives[n_mic+32] = C_gap_ferm;
  derivatives[n_mic+33] = C_gap_xa_dead;
  derivatives[n_mic+34] = C_gap_uptake;
  derivatives[n_mic+35] = C_gap_respir;
  derivatives[n_mic+36] = C_gap_urea;
  derivatives[n_mic+37] = N_gap_ferm;
  derivatives[n_mic+38] = N_gap_respir;

  // Other output, as a single flat vector (much faster for deSolve to unlist() than a nested list)
  // Names are created once in R (see outNames.R) and must match the order here
  int n_out = 9 + n_cf + 4 + 11 * n_mic + 1;
  NumericVector out(n_out);
  int k = 0;
  out[k++] = derivatives[n_mic + 23];  // COD_load_rate
  out[k++] = derivatives[n_mic + 24];  // C_load_rate
  out[k++] = derivatives[n_mic + 25];  // N_load_rate
  out[k++] = derivatives[n_mic + 17];  // CH4_emis_rate
  out[k++] = derivatives[n_mic + 18];  // CO2_emis_rate
  out[k++] = H2S_emis_rate;
  out[k++] = NH3_emis_rate_pit;
  out[k++] = NH3_emis_rate_floor;
  out[k++] = N2O_emis_rate;
  for (int i = 0; i < n_cf; ++i) out[k++] = interpTab(cf_tab, time_sum, i);
  for (int i = 0; i < n_mic; ++i) out[k++] = xa_fresh[i] * scale_xa_fresh;
  out[k++] = slurry_prod_rate;
  out[k++] = rain;
  out[k++] = evap;
  out[k++] = sum_rut;
  for (int i = 0; i < n_mic; ++i) out[k++] = rut[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = xa[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = qhat[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = qhat_opt[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = cum_inhib[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = HAC_inhib[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = NH3_inhib[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = NH4_inhib[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = H2S_inhib[i];
  for (int i = 0; i < n_mic; ++i) out[k++] = VFA / slurry_mass / ks[i];  // sat_index
  out[k++] = (VFA_H2_ferm + alpha_xa_dead * xa_dead) / sum_rut;  // VFA_supply_VFA_demand

  SEXP out_names = P(64);
  if (Rf_length(out_names) == n_out) {
    out.attr("names") = out_names;
  }

  return List::create(derivatives, out);
}
