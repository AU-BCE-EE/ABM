#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_derivatives(
  NumericVector yield, NumericVector rut, double slurry_prod_rate, double decay_rate, 
  NumericVector scale, NumericVector xa_fresh, NumericVector conc_fresh, 
  NumericVector ferm, NumericVector alpha, NumericVector COD_conv, 
  double rain, double evap, double area, double respiration, double sub_resp,
  double rut_urea, double NH3_emis_rate_pit, double NH3_emis_rate_floor, 
  double N2O_emis_rate, NumericVector rutsr, double H2S_emis_rate, 
  double VS_CH4, NumericVector lnA, NumericVector E_CH4, double R, 
  double temp_K, IntegerVector i_meth, double CO2_ferm_meth_sr
) {
  // Initialize output vector
  // alpha is:0 xa_dead,1 starch,2 Cfat,3 Cps,4 CPf,5 RFd,6 VSd,7 urea
  // conc_fresh is: 0 sulfide, 1 urea,2 sulfate,3 TAN,4 starch, 5 VFA,6 xa_aer, 7 xa_bac, 8 xa_dead, 9 Cfat, 10 CPs, 11 CPf, 12 RFd, 13 iNDF, 14 VSd, 15 VSd_A, 16 VSnd_A, 17 ash
  
  NumericVector derivatives(30);
  
  // Compute derivatives
  derivatives[0] = scale["yield"] * yield * rut + scale["xa_fresh"] * xa_fresh[0] * slurry_prod_rate - decay_rate * xa_fresh[0];
  derivatives[1] = slurry_prod_rate + (rain - evap) * area;
  derivatives[2] = slurry_prod_rate * conc_fresh["xa_aer"] + ferm["xa_aer_rate"] - decay_rate * conc_fresh["xa_aer"];
  derivatives[3] = slurry_prod_rate * conc_fresh["xa_bac"] + ferm["xa_bac_rate"] - decay_rate * conc_fresh["xa_bac"];
  
  derivatives[4] = slurry_prod_rate * conc_fresh["xa_dead"] - alpha["xa_dead"] * conc_fresh["xa_dead"] + 
    sum(decay_rate * xa_fresh) + decay_rate * (conc_fresh["xa_bac"] + conc_fresh["xa_aer"]) * 
    (1 - COD_conv["frac_CP_xa"]);
  
  derivatives[5] = slurry_prod_rate * conc_fresh["RFd"] - alpha["RFd"] * conc_fresh["RFd"] - 
    respiration * conc_fresh["RFd"] / sub_resp;
  
  derivatives[6] = slurry_prod_rate * conc_fresh["iNDF"];
  derivatives[7] = slurry_prod_rate * conc_fresh["ash"];
  
  derivatives[8] = slurry_prod_rate * conc_fresh["VSd"] - alpha["VSd"] * conc_fresh["VSd"] - 
    respiration * conc_fresh["VSd"] / sub_resp;
  
  derivatives[9] = slurry_prod_rate * conc_fresh["starch"] - alpha["starch"] * conc_fresh["starch"] - 
    respiration * conc_fresh["starch"] / sub_resp;
  
  derivatives[10] = slurry_prod_rate * conc_fresh["CPs"] - alpha["CPs"] * conc_fresh["CPs"] - 
    respiration * conc_fresh["CPs"] / sub_resp + decay_rate * 
    (conc_fresh["xa_bac"] + conc_fresh["xa_aer"]) * COD_conv["frac_CP_xa"];
  
  derivatives[11] = slurry_prod_rate * conc_fresh["CPf"] - alpha["CPf"] * conc_fresh["CPf"] - 
    respiration * conc_fresh["CPf"] / sub_resp;
  
  derivatives[12] = slurry_prod_rate * conc_fresh["Cfat"] - alpha["Cfat"] * conc_fresh["Cfat"] - 
    respiration * conc_fresh["Cfat"] / sub_resp;
  
  derivatives[13] = alpha["xa_dead"] * conc_fresh["xa_dead"] + ferm["VFA_H2"] - sum(rut) + 
    slurry_prod_rate * conc_fresh["VFA"];
  
  derivatives[14] = slurry_prod_rate * conc_fresh["urea"] - rut_urea;
  
  derivatives[15] = slurry_prod_rate * conc_fresh["TAN"] + rut_urea + ferm["TAN_min_ferm"] + 
    ferm["TAN_min_resp"] - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate;
  
  derivatives[16] = slurry_prod_rate * conc_fresh["sulfate"] - sum(rutsr) / COD_conv["S"];
  
  derivatives[17] = slurry_prod_rate * conc_fresh["sulfide"] + sum(rutsr) / COD_conv["S"] - H2S_emis_rate;
  
  derivatives[18] = -conc_fresh["VSd_A"] * (exp(lnA["VSd_A"] - E_CH4["VSd_A"] / (R * temp_K)) * 24 / 1000 * VS_CH4) + 
    slurry_prod_rate * conc_fresh["VSd_A"];
  
  derivatives[19] = slurry_prod_rate * conc_fresh["VSnd_A"];
  
  derivatives[20] = conc_fresh["VSd_A"] * (exp(lnA["VSd_A"] - E_CH4["VSd_A"] / (R * temp_K)) * 24 / 1000);
  
  derivatives[21] = NH3_emis_rate_pit + NH3_emis_rate_floor;
  derivatives[22] = N2O_emis_rate;
  
  derivatives[23] = sum(rut[i_meth]) / COD_conv["CH4"];
  derivatives[24] = CO2_ferm_meth_sr + ferm["CO2_resp"] + rut_urea / COD_conv["CO2_ureo"];
  
  derivatives[25] = sum(rut[i_meth]) + respiration + sum(rutsr);
  derivatives[26] = sum(rut[i_meth]);
  derivatives[27] = respiration;
  derivatives[28] = sum(rutsr);
  
  derivatives[29] = slurry_prod_rate * sum(conc_fresh[IntegerVector::create("starch", "VFA", "xa_aer", "xa_bac", "xa_dead", "Cfat", "CPs", "CPf", "RFd", "iNDF", "VSd")]) + 
    slurry_prod_rate * sum(xa_fresh * scale["xa_fresh"]);
  
  return derivatives;
}
