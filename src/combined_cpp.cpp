#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List combined_cpp(double g_NH4, bool pH_inhib_overrule,  // Single TRUE/FALSE
                   double pH, double pH_floor,
                   double TAN, double VFA, double sulfide, double slurry_mass,
                   NumericVector pH_LL, NumericVector pH_UL,
                   NumericVector ki_NH3_min, NumericVector ki_NH3_max,
                   NumericVector ki_NH4_min, NumericVector ki_NH4_max,
                   NumericVector ki_HAC, NumericVector ki_H2S_slope, 
                   NumericVector ki_H2S_int, NumericVector ki_H2S_min,
                   NumericVector IC50_low, 
                   double temp_K, double temp_C, double temp_standard, 
                   double area, double floor_area,
                   double Cfat, double CPs, double CPf, double RFd, 
                   double starch, double VSd, bool resp, List kl,
                   NumericVector qhat, NumericVector i_meth, NumericVector i_sr, NumericVector xa,
                   NumericVector ks_coefficient, double scale_ks, double ks_SO4, double sulfate,
                   double urea, double alpha_urea, double km_urea)
{
  // rough approximation from Rotz et al. IFSM 2012., "markfoged" thesis, "Petersen et al. 2014", "bilds?e et al. not published", "Elzing & Aarnik 1998", and own measurements..
  // Hard-wired equilibrium constants
  double log_HAC = -4.8288 + 21.42/temp_K;
  double log_NH3 = - 0.09046 - 2729.31/temp_K;
  double log_H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K);
  
  // Compute the fractions
  double HAC_frac = 1 - (1 / (1 + std::pow(10, -log_HAC - pH)));
  double NH3_frac = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH));
  double NH3_frac_floor = 1 / (1 + std::pow(10, -log_NH3 + std::log10(g_NH4) - pH_floor));
  double H2S_frac = 1 - (1 / (1 + std::pow(10, -log_H2S - pH)));
  
  int n = qhat.size();  // Number of microbial groups
  
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

// on change respiration if it is TRUE
if (resp && (slurry_mass / area >= 1)) {
  double kl_oxygen = exp(0.6158816 + 0.09205127 * temp_C);
  respiration = kl_oxygen * area * (kH_oxygen * 0.208) * (sub_resp / slurry_mass) / 100;
}

// ks is temp dependent
NumericVector ks = ks_coefficient * (0.8157 * exp(-0.063 * temp_C)); 
NumericVector rut(n);

// VFA consumption rate by sulfate reducers (g/d) affected by inhibition terms
for (int i = 0; i < i_sr.size(); ++i) {
  int idx_sr = i_sr[i];  // Convert from 1-based to 0-based index
  rut[idx_sr] = ((qhat[idx_sr] * VFA / slurry_mass * xa[idx_sr] / slurry_mass / (scale_ks * ks[idx_sr] + VFA / slurry_mass)) * slurry_mass *
                   (sulfate / slurry_mass) / (ks_SO4 + sulfate / slurry_mass)) * cum_inhib[idx_sr];
}

// VFA consumption rate by methanogen groups (g/d) affected by inhibition terms
for(int i = 0; i < i_meth.size(); ++i) {
  int idx = i_meth[i];
  rut[idx] = ((qhat[idx] * VFA / (slurry_mass) * xa[idx] / (slurry_mass)) / (scale_ks * ks[idx] + VFA / (slurry_mass)) * 
                (slurry_mass)) * cum_inhib[idx];
}

NumericVector rutsr(i_sr.size());

for (int i = 0; i < i_sr.size(); ++i) {
  int idx_sr = i_sr[i];  // Convert from 1-based to 0-based index
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

return List::create(Named("pH_inhib") = pH_inhib,
                    Named("NH3_inhib") = NH3_inhib,
                    Named("NH4_inhib") = NH4_inhib,
                    Named("HAC_inhib") = HAC_inhib,
                    Named("H2S_inhib") = H2S_inhib,
                    Named("cum_inhib") = cum_inhib,
                    Named("NH3_emis_rate_floor") = NH3_emis_rate_floor,
                    Named("NH3_emis_rate_pit") = NH3_emis_rate_pit,
                    Named("H2S_emis_rate") = H2S_emis_rate,
                    Named("respiration") = respiration,
                    Named("sub_resp") = sub_resp,
                    Named("rut_urea") = rut_urea,
                    Named("rut") = rut,
                    Named("rutsr") = rutsr);
}