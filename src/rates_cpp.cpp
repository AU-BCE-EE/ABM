#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
List rates_cpp(double t, NumericVector y, List parms, NumericVector p_idx, Rcpp::Function temp_C_fun, Rcpp::Function pH_fun, 
                  Rcpp::Function SO4_inhibition_fun, List conc_fresh_fun, NumericVector xa_fresh_fun, Rcpp::Function CTM_cpp){
  
  //y[y < 1E-10] <- 1E-10
  std::transform(y.begin(), y.end(), y.begin(), [](double value) {
    return (value < 1E-10) ? 1E-10 : value;
  });

  double slurry_prod_rate = parms[p_idx[0]];
  double t_run = parms[p_idx[1]];
  int n_mic = parms[p_idx[2]];
  
  //# NOT IMPLEMENTED YET
  //#if(!is.null(graze_int) & any(graze_int != 0)){
  //#  slurry_prod_rate <- graze_fun(t,  t_run, days, slurry_prod_rate, graze_int, graze_hours = graze[['hours_day']])
  //#}
  
  // not implemented variable conc_fresh and variable xa_fresh either
  // variable conc_frehs should be calculated from the List conc_fresh_fun
  // variable xa_fresh should be calculated from the xa_fresh_fun, which currently is a numeric vector or a data.frame.
  // need to always be of the same type (this could be a list or just an interpolation function).This is for later. 
  // print()Currently I pass it as a NumericVector (so variable is not an option currently)
  
  double pH = as<double>(pH_fun(t + t_run)); // missing the H2SO4 calc stuff.
  double temp_C = as<double>(temp_C_fun(t + t_run));
  double temp_K = temp_C + 273.15;
  double temp_standard = parms[p_idx[3]];
  
  // Extract state variable values from y argument
  NumericVector xa = y[Range(0, n_mic - 1)];
  
  double slurry_mass = y[n_mic];
  double xa_aer = y[n_mic+1];
  double xa_bac = y[n_mic+2];
  double xa_dead = y[n_mic+3];
  double RFd = y[n_mic+4];
  double iNDF = y[n_mic+5];
  double ash = y[n_mic+6];
  double VSd = y[n_mic+7];
  double starch = y[n_mic+8];
  double CPs = y[n_mic+9];
  double CPf = y[n_mic+10];
  double Cfat = y[n_mic+11];
  double VFA = y[n_mic+12];
  double urea = y[n_mic+13];
  double TAN = y[n_mic+14];
  double sulfate = y[n_mic+15];
  double sulfide = y[n_mic+16];
  double NH3_emis_cum = y[n_mic+17];
  double N2O_emis_cum = y[n_mic+18];
  double CH4_emis_cum = y[n_mic+19];
  double CO2_emis_cum = y[n_mic+20];
  double COD_conv_cum = y[n_mic+21];
  double COD_conv_cum_meth = y[n_mic+22];
  double COD_conv_cum_respir = y[n_mic+23];
  double COD_conv_cum_sr = y[n_mic+24];
  double COD_load_cum = y[n_mic+25];
  double C_load_cum = y[n_mic+26];
  double N_load_cum = y[n_mic+27];
  double slurry_load_cum = y[n_mic+28];

  // Arrhenius parms for hydrolysis stuff
  NumericVector A = parms[p_idx[4]];
  NumericVector E = parms[p_idx[5]];
  NumericVector scale = as<NumericVector>(parms[p_idx[6]]);
  
  double scale_ks_coefficient = scale[0];
  double scale_qhat_opt = scale[1];
  double scale_xa_fresh = scale[2];        
  double scale_yield = scale[3];
  double scale_alpha_opt = scale[4];
  
  double R = parms[p_idx[7]];
  double alpha_opt_scale_type = parms[p_idx[8]];
  double alpha_opt_scale_CP = parms[p_idx[9]];
  
  //# Hydrolysis rate with Arrhenius function cpp. 
  NumericVector alpha = A * exp(-E / (R * temp_K));
  alpha[Rcpp::Range(0,6)] = alpha[Rcpp::Range(0,6)] * scale_alpha_opt * alpha_opt_scale_type;
  alpha[Rcpp::Range(3,4)] = alpha[Rcpp::Range(3,4)] * alpha_opt_scale_CP;

  NumericVector T_opt = parms[p_idx[10]];
  NumericVector T_min = parms[p_idx[11]];
  NumericVector T_max = parms[p_idx[12]];
  NumericVector qhat_opt = parms[p_idx[13]];
  
  //# Microbial substrate utilization rate (vectorized calculation)
  NumericVector qhat = scale_qhat_opt * as<NumericVector>(CTM_cpp(temp_K, T_opt, T_min, T_max, qhat_opt));
  
  //# Decay of all microorganisms follow same kinetics with faster decay at higher temp up until 313, at which constant decay of 0.02
  NumericVector decay_rate_vec = parms[p_idx[14]];
  decay_rate_vec = as<NumericVector>(CTM_cpp(temp_K, 313, 273, 325, decay_rate_vec));
  double decay_rate_double = decay_rate_vec[0];
  
  if(temp_K > 313){
    decay_rate_double = 0.02;
  }
  
  double g_NH4 = parms[p_idx[15]];
  bool pH_inhib_overrule = parms[p_idx[16]];
  
  double pH_floor = parms[p_idx[17]];
  NumericVector pH_LL = parms[p_idx[18]];
  NumericVector pH_UL = parms[p_idx[19]];
  NumericVector ki_NH3_min = parms[p_idx[20]];
  NumericVector ki_NH3_max = parms[p_idx[21]];
  NumericVector ki_NH4_min = parms[p_idx[22]];
  NumericVector ki_NH4_max = parms[p_idx[23]];
  NumericVector ki_HAC = parms[p_idx[24]];
  NumericVector ki_H2S_slope = parms[p_idx[25]];
  NumericVector ki_H2S_int = parms[p_idx[26]];
  NumericVector ki_H2S_min = parms[p_idx[27]];
  NumericVector IC50_low = parms[p_idx[28]];
  double area = parms[p_idx[29]];
  double floor_area = parms[p_idx[30]];
  bool resp = parms[p_idx[31]];
  NumericVector kl = parms[p_idx[32]];
  NumericVector i_meth = parms[p_idx[33]];
  NumericVector i_sr = parms[p_idx[34]];
  NumericVector ks_coefficient = parms[p_idx[35]];
  double ks_SO4 = parms[p_idx[36]];
  
  double alpha_xa_dead = alpha[0];  // Assuming "xa_dead" is the first element
  double alpha_starch = alpha[1];   // Assuming "starch" is the second element
  double alpha_Cfat = alpha[2];     // And so on...
  double alpha_CPs = alpha[3];
  double alpha_CPf = alpha[4];
  double alpha_RFd = alpha[5];
  double alpha_VSd = alpha[6];
  double alpha_urea = alpha[7];
  
  
  double km_urea = parms[p_idx[37]];
  
  double log_HAC = -4.8288 + 21.42/temp_K;
  double log_NH3 = - 0.09046 - 2729.31/temp_K;
  double log_H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K);
  
  // Compute the fractions
  double HAC_frac = 1 - (1 / (1 + std::pow(10, -log_HAC - pH)));
  double NH3_frac = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH));
  double NH3_frac_floor = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH_floor));
  double H2S_frac = 1 - (1 / (1 + std::pow(10, -log_H2S - pH)));
  
  // use n_mic instead. int n = qhat.size();  // Number of microbial groups
  // inhib vector
  NumericVector pH_inhib(n_mic), NH3_inhib(n_mic), NH4_inhib(n_mic), HAC_inhib(n_mic), H2S_inhib(n_mic), cum_inhib(n_mic);
  
  if (pH_inhib_overrule) {
    // Compute pH inhibition and set cumulative inhibition = pH_inhib
    for (int i = 0; i < n_mic; ++i) {
      pH_inhib[i] = (1 + 2 * pow(10, 0.5 * (pH_LL[i] - pH_UL[i]))) /
        (1 + pow(10, pH - pH_UL[i]) + pow(10, pH_LL[i] - pH));
      cum_inhib[i] = pH_inhib[i];  // Directly set cumulative inhibition
    }
  } else {
    // Precompute invariant values (these do not change in the loop)
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
      
      // H2S inhibitionki_H2S_min
      double IC50 = (pH >= 6.8) ? (ki_H2S_slope[i] * pH + ki_H2S_int[i]) : IC50_low[i];
      double a = -0.5 / (IC50 - (H2S_frac * ki_H2S_min[i]));
      double b = 1 - (-0.5 / (IC50 - (H2S_frac * ki_H2S_min[i])) * H2S_frac * ki_H2S_min[i] / slurry_mass);
      
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
  
  // Extract values from kl list
  double kl_NH3 = kl[0];
  double kl_NH3_floor = kl[1];
  double kl_H2S = kl[2];
  
  // NH3 emission rates
  double NH3_emis_rate_floor = kl_NH3_floor * floor_area * ((NH3_frac_floor * TAN) / slurry_mass) * 1000 / H_NH3;
  double NH3_emis_rate_pit = kl_NH3 * area * ((NH3_frac * TAN) / slurry_mass) * 1000 / H_NH3;
  
  // H2S emission rate
  double H2S_emis_rate = kl_H2S * area * ((H2S_frac * sulfide) / slurry_mass) * 1000;
  
  // Respiration calculation
  double respiration = 0.0;
  double sub_resp = Cfat + CPs + CPf + RFd + starch + VSd;
  
  // on change respiration if it is TRUE
  if (resp && (slurry_mass / area >= 1)) {
    double kl_oxygen = exp(0.6158816 + 0.09205127 * temp_C);
    respiration = kl_oxygen * area * (kH_oxygen * 0.208) * (sub_resp / slurry_mass) / 100;
  }
  
  // ks is temp dependent
  NumericVector ks = ks_coefficient * (0.8157 * exp(-0.063 * temp_C)); 
  NumericVector rut(n_mic);
  
  // VFA consumption rate by sulfate reducers (g/d) affected by inhibition terms
  for (int i = 0; i < i_sr.size(); ++i) {
    int idx_sr = i_sr[i];  // Convert from 1-based to 0-based index
    rut[idx_sr] = ((qhat[idx_sr] * VFA / slurry_mass * xa[idx_sr] / slurry_mass / (scale_ks_coefficient * ks[idx_sr] + VFA / slurry_mass)) * slurry_mass *
      (sulfate / slurry_mass) / (ks_SO4 + sulfate / slurry_mass)) * cum_inhib[idx_sr];
  }
  
  // VFA consumption rate by methanogen groups (g/d) affected by inhibition terms
  for(int i = 0; i < i_meth.size(); ++i) {
    int idx = i_meth[i];
    rut[idx] = ((qhat[idx] * VFA / (slurry_mass) * xa[idx] / (slurry_mass)) / (scale_ks_coefficient * ks[idx] + VFA / (slurry_mass)) * 
      (slurry_mass)) * cum_inhib[idx];
  }
  
  NumericVector rutsr(i_sr.size());
  
  for (int i = 0; i < i_sr.size(); ++i) {
    int idx_sr = i_sr[i];  // Convert from 1-based to 0-based index
    if (idx_sr >= 0 && idx_sr < n_mic) {
      rutsr[i] = rut[idx_sr];
    }
  }
  
  // If rutsr is empty, initialize with a value of 0
  if (rutsr.size() == 0) {
    double rutsr = 0;  // Initialize with 0 if empty
  }
  
  // stop if any rut are below 0
  if (std::any_of(rut.begin(), rut.end(), [](double val) { return val < 0; })) {
    Rcpp::stop("In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)");
  }
  // urea hydrolysis by Michaelis Menten
  double rut_urea = (alpha_urea * urea) / (km_urea + urea/slurry_mass);
  
  //# CO2 production from fermentation + methanogenesis + sulfate reduction + aerobic respiration at slurry surface
  //# also calcualtes growth rate of xa_bac and xa_aer, mineralization rates and COD production from fermentation
  
  //### MOVE all the stuff below until derivatives to CPP?
  //ferm <- stoich(alpha, y, conc_fresh, sub_resp, respiration,
  //               carb, pro, lip, carb_resp, pro_resp, lip_resp, ace, hyd,
  //               ace_sr, hyd_sr)
  
  //CO2_ferm_meth_sr <- ferm$ferm[['CO2']] * 44.01 + sum(rut[i_meth+1])/ferm$COD_conv_meth_CO2[[1]] + sum(rutsr)/ferm$COD_conv_sr_CO2[[1]]
  //CO2_ferm <- ferm$ferm[['CO2']] * 44.01 
  
  
  // derivatives
  
  List conc_fresh = parms[p_idx[38]];
  
  double conc_fresh_sulfide = conc_fresh[0];  // Assuming "sulfide" is the first element
  double conc_fresh_urea = conc_fresh[1];     // Assuming "urea" is the second element
  double conc_fresh_sulfate = conc_fresh[2];  // And so on...
  double conc_fresh_TAN = conc_fresh[3];
  double conc_fresh_starch = conc_fresh[4];
  double conc_fresh_VFA = conc_fresh[5];
  double conc_fresh_xa_aer = conc_fresh[6];
  double conc_fresh_xa_bac = conc_fresh[7];
  double conc_fresh_xa_dead = conc_fresh[8];
  double conc_fresh_Cfat = conc_fresh[9];
  double conc_fresh_CPs = conc_fresh[10];
  double conc_fresh_CPf = conc_fresh[11];
  double conc_fresh_RFd = conc_fresh[12];
  double conc_fresh_iNDF = conc_fresh[13];
  double conc_fresh_VSd = conc_fresh[14];
  double conc_fresh_ash = conc_fresh[15];
  
  // NOT IMPLEMENTED YET. End products of fermentation (previously from stoich()). 
  // Right now are just placeholder to test speed. Need to be properly calculated in hard_pars() or otherwise
  double xa_bac_rate = parms[p_idx[39]];
  double TAN_min_ferm = parms[p_idx[40]];
  double VFA_H2_ferm = parms[p_idx[41]];
  double COD_conv_meth_CO2 = parms[p_idx[42]];
  double COD_conv_sr_CO2 = parms[p_idx[43]];
  double CO2_ferm = parms[p_idx[44]];
  
  // below are calculated from conc_fresh composition  
  double TAN_min_resp = as<double>(parms[p_idx[45]]) * respiration;
  double CO2_resp = as<double>(parms[p_idx[46]]) * respiration;
  double xa_aer_rate = as<double>(parms[p_idx[47]]) * respiration;
  
  NumericVector COD_conv = parms[p_idx[48]];
  
  double COD_conv_frac_CP_xa = COD_conv[24];
  double COD_conv_S = COD_conv[8];
  double COD_conv_CH4 = COD_conv[0];
  double COD_conv_CO2_ureo = COD_conv[12];
  double COD_conv_CP_N = COD_conv[13];
  
  double COD_conv_C_starch = COD_conv[18];
  double COD_conv_C_VFA = COD_conv[21];
  double COD_conv_C_xa = COD_conv[14];
  double COD_conv_C_Cfat = COD_conv[19];
  double COD_conv_C_CP = COD_conv[20];
  double COD_conv_C_VSd = COD_conv[22];
  double COD_conv_C_N_urea = COD_conv[23];
  
  double sum_rut = std::accumulate(rut.begin(), rut.end(), 0.0);
  double sum_rutmeth = 0.0;
  for(int i = 0; i < i_meth.size(); ++i) {
    int idx = i_meth[i];
    sum_rutmeth += rut[idx];
  }
  
  double sum_rutsr = 0.0;
  for(int i = 0; i < i_sr.size(); ++i) {
    int idxsr = i_sr[i];
    sum_rutsr += rut[idxsr];
  }
  
  // this calculates the CO2 produced from combined fermentation, methanogenesis and sulfate reduction
  double CO2_ferm_meth_sr = CO2_ferm * 44.01 + sum_rutmeth/COD_conv_meth_CO2 + sum_rutsr/COD_conv_sr_CO2;
  
  NumericVector xa_fresh = parms[p_idx[49]];
  
  double sum_xa_fresh = std::accumulate(xa_fresh.begin(), xa_fresh.end(), 0.0);

  int nn = n_mic + 29;
  
  NumericVector derivatives(nn);
  
  NumericVector yield = parms[p_idx[50]];
  
  for (int i = 0; i < n_mic; i++) {
    derivatives[i] = scale_yield * yield[i] * rut[i] + scale_xa_fresh * xa_fresh[i] * slurry_prod_rate - decay_rate_double * xa[i];
  }
  
  double rain = parms[p_idx[51]];
  double evap = parms[52];
  double N2O_emis_rate = parms[p_idx[53]];
  
  derivatives[n_mic] = slurry_prod_rate + (rain - evap) * area;
  derivatives[n_mic+1] = slurry_prod_rate * conc_fresh_xa_aer + xa_aer_rate - decay_rate_double * xa_aer;
  derivatives[n_mic+2] = slurry_prod_rate * conc_fresh_xa_bac + xa_bac_rate - decay_rate_double * xa_bac;
  derivatives[n_mic+3] = slurry_prod_rate * conc_fresh_xa_dead - alpha_xa_dead * xa_dead + decay_rate_double * std::accumulate(xa.begin(), xa.end(), 0.0) + decay_rate_double * (xa_bac + xa_aer) * (1 - COD_conv_frac_CP_xa);
  derivatives[n_mic+4] = slurry_prod_rate * conc_fresh_RFd - alpha_RFd * RFd - respiration * RFd / sub_resp;
  derivatives[n_mic+5] = slurry_prod_rate * conc_fresh_iNDF;
  derivatives[n_mic+6] = slurry_prod_rate * conc_fresh_ash;
  derivatives[n_mic+7] = slurry_prod_rate * conc_fresh_VSd - alpha_VSd * VSd - respiration * VSd / sub_resp;
  derivatives[n_mic+8] = slurry_prod_rate * conc_fresh_starch - alpha_starch * starch - respiration * starch / sub_resp;
  derivatives[n_mic+9] = slurry_prod_rate * conc_fresh_CPs - alpha_CPs * CPs - respiration * CPs / sub_resp + decay_rate_double * (xa_bac + xa_aer) * COD_conv_frac_CP_xa;
  derivatives[n_mic+10] = slurry_prod_rate * conc_fresh_CPf - alpha_CPf * CPf - respiration * CPf / sub_resp;
  derivatives[n_mic+11] = slurry_prod_rate * conc_fresh_Cfat - alpha_Cfat * Cfat - respiration * Cfat / sub_resp;
  derivatives[n_mic+12] = alpha_xa_dead * xa_dead + VFA_H2_ferm - sum_rut + slurry_prod_rate * conc_fresh_VFA;
  derivatives[n_mic+13] = slurry_prod_rate * conc_fresh_urea - rut_urea;
  derivatives[n_mic+14] = slurry_prod_rate * conc_fresh_TAN + rut_urea + TAN_min_ferm + TAN_min_resp - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate;
  derivatives[n_mic+15] = slurry_prod_rate * conc_fresh_sulfate - sum_rutsr / COD_conv_S;
  derivatives[n_mic+16] = slurry_prod_rate * conc_fresh_sulfide + sum_rutsr / COD_conv_S - H2S_emis_rate;
  derivatives[n_mic+17] = NH3_emis_rate_pit + NH3_emis_rate_floor;
  derivatives[n_mic+18] = N2O_emis_rate;
  derivatives[n_mic+19] = sum_rutmeth / COD_conv_CH4;
  derivatives[n_mic+20] = CO2_ferm_meth_sr + CO2_resp + rut_urea / COD_conv_CO2_ureo;
  derivatives[n_mic+21] = sum_rut + respiration;
  derivatives[n_mic+22] = sum_rutmeth;
  derivatives[n_mic+23] = respiration;
  derivatives[n_mic+24] = sum_rutsr;
  derivatives[n_mic+25] = slurry_prod_rate * (conc_fresh_starch + conc_fresh_VFA + conc_fresh_xa_aer + conc_fresh_xa_bac + conc_fresh_xa_dead + conc_fresh_Cfat + conc_fresh_CPs + conc_fresh_CPf + conc_fresh_RFd + conc_fresh_iNDF + conc_fresh_VSd + (sum_xa_fresh * scale_xa_fresh));
  derivatives[n_mic+26] = slurry_prod_rate * ((conc_fresh_starch + conc_fresh_RFd + conc_fresh_iNDF)/COD_conv_C_starch + conc_fresh_VFA/COD_conv_C_VFA + (conc_fresh_xa_aer + conc_fresh_xa_bac + conc_fresh_xa_dead)/COD_conv_C_xa + 
    conc_fresh_Cfat/COD_conv_C_Cfat + (conc_fresh_CPs + conc_fresh_CPf)/COD_conv_C_CP + 
    conc_fresh_VSd/COD_conv_C_VSd + conc_fresh_urea/COD_conv_C_N_urea + (sum_xa_fresh * scale_xa_fresh)/(COD_conv_C_xa));
  derivatives[n_mic+27] =  slurry_prod_rate * ((conc_fresh_CPs + conc_fresh_CPf)/COD_conv_CP_N + conc_fresh_TAN + conc_fresh_urea);
  derivatives[n_mic+28] = slurry_prod_rate;
 
    // Return the main list containing both derivatives and rates
    return List::create(Named("derivatives") = derivatives,
                        Named("COD_load_rate") = derivatives[n_mic + 25],
                        Named("C_load_rate") = derivatives[n_mic + 26],
                        Named("N_load_rate") = derivatives[n_mic + 27],
                        Named("CH4_emis_rate") = derivatives[n_mic + 19],
                        Named("CO2_emis_rate") = derivatives[n_mic + 20],
                        Named("H2S_emis_rate") = H2S_emis_rate,
                        Named("NH3_emis_rate_pit") = NH3_emis_rate_pit,
                        Named("NH3_emis_rate_floor") = NH3_emis_rate_floor,
                        Named("N2O_emis_rate") = N2O_emis_rate,
                        Named("qhat") = qhat,
                        Named("alpha") = alpha,
                        Named("H2S_inhib") = H2S_inhib,
                        Named("NH3_inhib") = NH3_inhib,
                        Named("NH4_inhib") = NH4_inhib,
                        Named("HAC_inhib") = HAC_inhib,
                        Named("cum_inhib") = cum_inhib,
                        Named("conc_fresh") = conc_fresh,
                        Named("xa_fresh") = xa_fresh * scale_xa_fresh,
                        Named("area") = area,
                        Named("slurry_prod_rate") = slurry_prod_rate,
                        Named("respiration") = respiration,  
                        Named("rain") = rain,
                        Named("evap") = evap);
  }
  

