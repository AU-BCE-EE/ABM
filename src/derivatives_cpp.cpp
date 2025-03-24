#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector derivatives_cpp(
  NumericVector yield, NumericVector rut, double slurry_prod_rate, double decay_rate, 
  NumericVector scale, NumericVector xa_fresh, List conc_fresh, 
  List ferm, NumericVector alpha, NumericVector COD_conv, 
  double rain, double evap, double area, double respiration, double sub_resp,
  double rut_urea, double NH3_emis_rate_pit, double NH3_emis_rate_floor, 
  double N2O_emis_rate, NumericVector rutsr, double H2S_emis_rate, 
  double R, double temp_K, NumericVector i_meth, NumericVector i_sr, double CO2_ferm_meth_sr, 
  NumericVector xa, double slurry_mass, double xa_bac, double xa_aer, double xa_dead,
  double RFd, double iNDF, double ash, double VSd, double starch, double CPs, double CPf,
  double Cfat, double VFA, double urea, double TAN, double sulfate, double sulfide, NumericVector qhat
) {
  
  // Initialize output vector
  double alpha_xa_dead = alpha[0];  // Assuming "xa_dead" is the first element
  double alpha_starch = alpha[1];   // Assuming "starch" is the second element
  double alpha_Cfat = alpha[2];     // And so on...
  double alpha_CPs = alpha[3];
  double alpha_CPf = alpha[4];
  double alpha_RFd = alpha[5];
  double alpha_VSd = alpha[6];

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

  double scale_xa_fresh = scale[2];        // And so on...
  double scale_yield = scale[3];
  
  double ferm_xa_bac_rate = ferm[3];
  double ferm_xa_aer_rate = ferm[4];
  double ferm_TAN_min_ferm = ferm[5];
  double ferm_TAN_min_resp = ferm[6];
  double ferm_VFA_H2 = ferm[7];
  double ferm_CO2_resp = ferm[8];
  
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

  double sum_xa_fresh = std::accumulate(xa_fresh.begin(), xa_fresh.end(), 0.0);
  double n = qhat.size();
  int nn = n + 29;
  
  NumericVector derivatives(nn);
  
  for (int i = 0; i < n; i++) {
    derivatives[i] = scale_yield * yield[i] * rut[i] + scale_xa_fresh * xa_fresh[i] * slurry_prod_rate - decay_rate * xa[i];
  }

  derivatives[n] = slurry_prod_rate + (rain - evap) * area;
  derivatives[n+1] = slurry_prod_rate * conc_fresh_xa_aer + ferm_xa_aer_rate - decay_rate * xa_aer;
  derivatives[n+2] = slurry_prod_rate * conc_fresh_xa_bac + ferm_xa_bac_rate - decay_rate * xa_bac;
  derivatives[n+3] = slurry_prod_rate * conc_fresh_xa_dead - alpha_xa_dead * xa_dead + decay_rate * std::accumulate(xa.begin(), xa.end(), 0.0) + decay_rate * (xa_bac + xa_aer) * (1 - COD_conv_frac_CP_xa);
  derivatives[n+4] = slurry_prod_rate * conc_fresh_RFd - alpha_RFd * RFd - respiration * RFd / sub_resp;
  derivatives[n+5] = slurry_prod_rate * conc_fresh_iNDF;
  derivatives[n+6] = slurry_prod_rate * conc_fresh_ash;
  derivatives[n+7] = slurry_prod_rate * conc_fresh_VSd - alpha_VSd * VSd - respiration * VSd / sub_resp;
  derivatives[n+8] = slurry_prod_rate * conc_fresh_starch - alpha_starch * starch - respiration * starch / sub_resp;
  derivatives[n+9] = slurry_prod_rate * conc_fresh_CPs - alpha_CPs * CPs - respiration * CPs / sub_resp + decay_rate * (xa_bac + xa_aer) * COD_conv_frac_CP_xa;
  derivatives[n+10] = slurry_prod_rate * conc_fresh_CPf - alpha_CPf * CPf - respiration * CPf / sub_resp;
  derivatives[n+11] = slurry_prod_rate * conc_fresh_Cfat - alpha_Cfat * Cfat - respiration * Cfat / sub_resp;
  derivatives[n+12] = alpha_xa_dead * xa_dead + ferm_VFA_H2 - sum_rut + slurry_prod_rate * conc_fresh_VFA;
  derivatives[n+13] = slurry_prod_rate * conc_fresh_urea - rut_urea;
  derivatives[n+14] = slurry_prod_rate * conc_fresh_TAN + rut_urea + ferm_TAN_min_ferm + ferm_TAN_min_resp - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate;
  derivatives[n+15] = slurry_prod_rate * conc_fresh_sulfate - sum_rutsr / COD_conv_S;
  derivatives[n+16] = slurry_prod_rate * conc_fresh_sulfide + sum_rutsr / COD_conv_S - H2S_emis_rate;
  derivatives[n+17] = NH3_emis_rate_pit + NH3_emis_rate_floor;
  derivatives[n+18] = N2O_emis_rate;
  derivatives[n+19] = sum_rutmeth / COD_conv_CH4;
  derivatives[n+20] = CO2_ferm_meth_sr + ferm_CO2_resp + rut_urea / COD_conv_CO2_ureo;
  derivatives[n+21] = sum_rut + respiration;
  derivatives[n+22] = sum_rutmeth;
  derivatives[n+23] = respiration;
  derivatives[n+24] = sum_rutsr;
  derivatives[n+25] = slurry_prod_rate * (conc_fresh_starch + conc_fresh_VFA + conc_fresh_xa_aer + conc_fresh_xa_bac + conc_fresh_xa_dead + conc_fresh_Cfat + conc_fresh_CPs + conc_fresh_CPf + conc_fresh_RFd + conc_fresh_iNDF + conc_fresh_VSd + (sum_xa_fresh * scale_xa_fresh));
  derivatives[n+26] = slurry_prod_rate * ((conc_fresh_starch + conc_fresh_RFd + conc_fresh_iNDF)/COD_conv_C_starch + conc_fresh_VFA/COD_conv_C_VFA + (conc_fresh_xa_aer + conc_fresh_xa_bac + conc_fresh_xa_dead)/COD_conv_C_xa + 
    conc_fresh_Cfat/COD_conv_C_Cfat + (conc_fresh_CPs + conc_fresh_CPf)/COD_conv_C_CP + 
    conc_fresh_VSd/COD_conv_C_VSd + conc_fresh_urea/COD_conv_C_N_urea + (sum_xa_fresh * scale_xa_fresh)/(COD_conv_C_xa));
  derivatives[n+27] =  slurry_prod_rate * ((conc_fresh_CPs + conc_fresh_CPf)/COD_conv_CP_N + conc_fresh_TAN + conc_fresh_urea);
  derivatives[n+28] = slurry_prod_rate;

  return derivatives;
}
