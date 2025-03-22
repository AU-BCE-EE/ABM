#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List inhib_cpp(
    bool pH_inhib_overrule,  // Single TRUE/FALSE
    double pH, double NH3_frac, double HAC_frac, double H2S_frac,
    double TAN, double VFA, double sulfide, double slurry_mass,
    NumericVector pH_LL, NumericVector pH_UL,
    NumericVector ki_NH3_min, NumericVector ki_NH3_max,
    NumericVector ki_NH4_min, NumericVector ki_NH4_max,
    NumericVector ki_HAC, NumericVector ki_H2S_slope, 
    NumericVector ki_H2S_int, NumericVector ki_H2S_min,
    NumericVector IC50_low
) {
  int n = pH_LL.size();  // Number of microbial groups
  
  // Output vectors
  NumericVector pH_inhib(n), NH3_inhib(n), NH4_inhib(n), HAC_inhib(n), H2S_inhib(n), cum_inhib(n);
  
  if (pH_inhib_overrule) {
    // Compute pH inhibition and set cumulative inhibition = pH_inhib
    for (int i = 0; i < n; ++i) {
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
    
    for (int i = 0; i < n; ++i) {
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
  
  return List::create(
    _["pH_inhib"] = pH_inhib,
    _["NH3_inhib"] = NH3_inhib,
    _["NH4_inhib"] = NH4_inhib,
    _["HAC_inhib"] = HAC_inhib,
    _["H2S_inhib"] = H2S_inhib,
    _["cum_inhib"] = cum_inhib
  );
}
