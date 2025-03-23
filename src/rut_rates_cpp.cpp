#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rut_rates_cpp(double temp_K, double temp_C, double temp_standard, double slurry_mass, 
                    double area, double floor_area, double NH3_frac, double NH3_frac_floor, double TAN,
                    double H2S_frac, double sulfide, double Cfat, double CPs, double CPf, double RFd, 
                    double starch, double VSd, bool resp, List kl,
                    NumericVector qhat, NumericVector i_meth, NumericVector i_sr, NumericVector xa,
                    double VFA, double scale_ks, NumericVector ks, double ks_SO4, double sulfate,
                    NumericVector cum_inhib, double urea, double alpha_urea, double km_urea)
  {
  
  // Henry's constant for NH3
  double H_NH3 = 1431 * pow(1.053, (293 - temp_K));
  
  // Henry's constant for oxygen
  double kH_oxygen = 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000;
  
  // Extract values from kl list
  double kl_NH3_floor = as<double>(kl["NH3_floor"]);
  double kl_NH3 = as<double>(kl["NH3"]);
  double kl_H2S = as<double>(kl["H2S"]);
  
  // NH3 emission rates
  double NH3_emis_rate_floor = kl_NH3_floor * floor_area * ((NH3_frac_floor * TAN) / slurry_mass) * 1000 / H_NH3;
  double NH3_emis_rate_pit = kl_NH3 * area * ((NH3_frac * TAN) / slurry_mass) * 1000 / H_NH3;
  
  // H2S emission rate
  double H2S_emis_rate = kl_H2S * area * ((H2S_frac * sulfide) / slurry_mass) * 1000;
  
  // Respiration calculation
  double respiration = 0.0;
  double sub_resp = Cfat + CPs + CPf + RFd + starch + VSd;
  
  if (resp && (slurry_mass / area >= 1)) {
    double kl_oxygen = exp(0.6158816 + 0.09205127 * temp_C);
    respiration = kl_oxygen * area * (kH_oxygen * 0.208) * (sub_resp / slurry_mass) / 100;
  }
  
  int n = qhat.size();
  
  NumericVector rut(n);
  // VFA consumption rate by sulfate reducers (g/d) affected by inhibition terms
  for (int i = 0; i < i_sr.size(); ++i) {
    int idx_sr = i_sr[i] - 1;  // Convert from 1-based to 0-based index
      rut[idx_sr] = ((qhat[idx_sr] * VFA / slurry_mass * xa[idx_sr] / slurry_mass / (scale_ks * ks[idx_sr] + VFA / slurry_mass)) * slurry_mass *
        (sulfate / slurry_mass) / (ks_SO4 + sulfate / slurry_mass)) * cum_inhib[idx_sr];
  }
  
  // VFA consumption rate by methanogen groups (g/d) affected by inhibition terms
  for(int i = 0; i < i_meth.size(); ++i) {
    int idx = i_meth[i] - 1;
      rut[idx] = ((qhat[idx] * VFA / (slurry_mass) * xa[idx] / (slurry_mass)) / (scale_ks * ks[idx] + VFA / (slurry_mass)) * 
        (slurry_mass)) * cum_inhib[idx];
  }
  
  NumericVector rutsr(i_sr.size());
    for (int i = 0; i < i_sr.size(); ++i) {
      int idx_sr = i_sr[i] - 1;  // Convert from 1-based to 0-based index
        if (idx_sr >= 0 && idx_sr < rut.size()) {
        rutsr[i] = rut[idx_sr];
        }
    }
  
  // If rutsr is empty, initialize with a value of 0
  if (rutsr.size() == 0) {
    double rutsr = 0;  // Initialize with 0 if empty
  }
  
  if (std::any_of(rut.begin(), rut.end(), [](double val) { return val < 0; })) {
    Rcpp::stop("In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)");
  }
  // urea hydrolysis by Michaelis Menten
  double rut_urea = (alpha_urea * urea) / (km_urea + urea/slurry_mass);
  
  return List::create(Named("NH3_emis_rate_floor") = NH3_emis_rate_floor,
                      Named("NH3_emis_rate_pit") = NH3_emis_rate_pit,
                      Named("H2S_emis_rate") = H2S_emis_rate,
                      Named("respiration") = respiration,
                      Named("sub_resp") = sub_resp,
                      Named("rut_urea") = rut_urea,
                      Named("rut") = rut,
                      Named("rutsr") = rutsr);
                      
}
