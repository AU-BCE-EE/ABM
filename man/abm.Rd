\name{abm}
\alias{abm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The ABM model
}
\description{
  An implementation of the ABM model for simulation of anaerobic biodegradation of organic substrates.
}
\usage{
abm(days = 365,
    delta_t = 1,
    mng_pars = list(slurry_prod_rate = 1000,
                    slurry_mass = 0,
                    max_slurry_mass = 33333,
                    resid_frac = 0.10,
                    area = 11,
                    temp_C = 20),
    man_pars = list(conc_fresh = list(S2 = 0.0, SO4 = 0.2, TAN = 1.0, VFA = 4.2, Sp = 65, COD = 160),
                    pH = 7),
   grp_pars = list(grps = c('m1', 'm2', 'm3', 'm4', 'm5'),
                   yield = c(all = 0.04),
                   xa_fresh = c(default = 0.001, m3 = 0.01),
                   xa_init = c(m1 = 0.01, m2 = 0.005, m3 = 0.005, m4 = 0.005, m5 = 0.001),
                   decay_rate = c(all = 0.02),
                   ks_coefficient = c(all = 1.0),
                   resid_enrich = c(all = 0.0),
                   qhat_opt = c(m1 = 3.6, m2 = 5.6 , m3 = 7.2, m4 = 8, m5 = 8),
                   T_opt = c(m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55),
                   T_min = c(m1 = 0, m2 = 8, m3 = 15, m4 = 26.25, m5 = 30),
                   T_max = c(m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60),
                   ki_NH3_min = c(all = 0.015),
                   ki_NH3_max = c(all = 0.13),
                   ki_NH4_min = c(all = 2.7),
                   ki_NH4_max = c(all = 4.8),
                   pH_upr = c(all = 8.0),
                   pH_lwr = c(all = 6.5)),
    mic_pars = list(ks_SO4 = 0.0067,
                    ki_H2S_meth = 0.23,
                    ki_H2S_sr = 0.25,
                    alpha_opt = 0.015,
                    alpha_T_opt = 313,
                    alpha_T_min = 273,
                    alpha_T_max = 320.67),
    chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.57, CO2_aer = 1.3, CO2_sr = 1.3), kl = c(H2S = 0.02, oxygen = 0.5)),
    add_pars = NULL,
    startup = -Inf,
    starting = NULL,
    approx_method_temp = 'linear',
    approx_method_pH = 'linear',
    approx_method_SO4 = 'linear',
    par_key = '\\.',
    value = 'ts',
    warn = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{days}{
    Total number of days in simulation.
}
  \item{delta_t}{
    Time step in output (d).
}
  \item{mng_pars}{
    List of management parameters including slurry production rate (kg/day) (\code{slurry_prod_rate}),
    initial slurry mass in the system (kg) (\code{slurry_mass}),
    maximum slurry capacity of the system (kg) (\code{max_slurry_mass}),
    fraction of slurry retained in the system after being emptied (\code{resid_frac}),
    surface area of the slurry in the system (m^2) (\code{area}),
    slurry temperature (deg. C) (\code{temp_C}).
    Each element may be a length-one vector or, in some cases, a list.
    See usage, examples, and vignettes.
}
  \item{man_pars}{
    List of manure parameters.
    Elements are:
    \code{conc_fresh}, the concentration of slurry components in fresh slurry;
    \code{pH}, pH of the slurry entering the system and in the slurry in the system.
    The \code{conc_fresh} element is a list itself, with the following elements:
    \code{S2} (sulfide), \code{SO4} (sulfate), \code{TAN} (ammonia nitrogen), \code{VFA} (volatile fatty acids), \code{Sp} (degradable, particulate, organic substrate), \code{COD} (total COD).
    Concentrations are on an elemental basis (N and S) or else COD (organic components).
}
  \item{grp_pars}{
   List containing parameter values for the microbial groups included in the model.
   Elements are:
   \code{yield}, Growth yield of microbial groups (g COD-biomass / g COD-substrate) ;
   \code{xa_fresh}, Concentration of microbial groups in the slurry entering the system (g COD-biomass/kg);
   \code{xa_init}, Initial concentration of microbial groups in the slurry (g COD-biomass/kg);
   \code{decay_rate}, Deacay rate of microbial groups (1/day);
   \code{ks_coefficient}, Coefficient used to tune the half saturation constant (g COD-substrate/kg);
   \code{resid_enrich}, Enrichment of microbial groups in the slurry after being emptied ;
   \code{qhat_opt}, Maximum specific substrate utilization rate of microbial groups at their optimum growth temperature (g COD-substrate/g COD-biomass/day);
   \code{T_opt}, Optimum growth temperature of microbial groups (K);
   \code{T_min}, Minimum viable growth temperature of microbial groups (K);
   \code{T_max}, Maximum vaible growth temperature of microbial groups (K);
   \code{ki_NH3_min}, NH3 concentration at which NH3 inhibition starts (g TAN-nitrogen/kg);
   \code{ki_NH3_max}, NH3 concentration wher microbial groups are 100 percent inhibited (g TAN-nitrogen/kg);
   \code{ki_NH4_min}, NH4 concentration at which NH3 inhibition starts (g TAN-nitrogen/kg);
   \code{ki_NH4_max}, NH4 concentration wher microbial groups are 100 percent inhibited (g TAN-nitrogen/kg);
   \code{pH_upr}, upper pH where 50 percent inhibition of microbial groups occurs;
   \code{pH_upr}, lower pH where 50 percent inhibition of microbial groups occurs.
   Each of these elements is a vector with a named element for each group.
   See usage, examples, and vignettes.
}
  \item{mic_pars}{
    List of other microbial parameters.
    Elements are:
    \code{ks_SO4}, Half sulfate saturation constant for sulfate reducers (g sulfate-sulfur/kg);
    \code{ki_H2S_meth}, hydrogen sulfide concentration at which 100 percent inhibition of methanogen groups occur (g sulfide-sulfur/kg);
    \code{ki_H2S_sr}, hydrogen sulfide concentration at which 100 percent inhibition of sulfate reducer groups occur (g sulfide-sulfur/kg);
    \code{alpha_opt}, hydrolysis rate constant at the optimum hydrolysis temperature (1/day);
    \code{alpha_T_opt}, optimum temperature for hydrolysis (K);
    \code{alpha_T_min}, minimum temperature where hydrolysis occurs (K);
    \code{alpha_T_max}, maximum temperature where hydrolysis occurs (K).
}
  \item{chem_pars}{
    List of other chemical parameters.
    Elements are:
    \code{kl}, Mass transfer coefficient of oxygen and hydrogen sulfide through the slurry/air interface, liquid-phase units (m/day);
    \code{COD_conv}, conversion factors for using a COD mass basis with the following named elements:
    \code{CH4}, Methane productivity coefficient (g CH4/g COD-substrate) ;
    \code{S}, COD-substrate to sulfur (g sulfur/g COD-substrate) ;
    \code{VS}, COD to volatile solids (VS) (g VS/g COD), ;
    \code{CO2_anaer}, Anaerobic carbon dioxide productivity coefficient (g CO2/g COD-substrate);
    \code{CO2_aer}, Aerobic carbon dioxide productivity coefficient (g CO2/g COD-substrate);
    \code{CO2_sr}, carbon dioxide productivity coefficient for sulfate reduction (g CO2/g COD-substrate).
}
  \item{add_pars}{
    List of other additional parameters not set in other arguments, or to be modified from defaults.
    See vignette.
}
  \item{startup}{
    A ``startup'' period that should be excluded from the results.
    Length-one numeric vector (d).
    Optional.
}
  \item{starting}{
    Starting conditions.
    Output (\code{value = 'ts'}) from a previous \code{abm()} call.
    Values here for microbial biomass and slurry mass override those set in other arguments.
    Optional.
}
  \item{approx_method_temp}{
    Sets type of interpolation when temperature is variable (\code{temp_c} is a data frame).
    Passed directly to \code{\link{approxfun}}.
    See vignette.
}
  \item{approx_method_pH}{
    Sets type of interpolation when pH is variable (\code{pH} is a data frame).
    Passed directly to \code{\link{approxfun}}.
    See vignette.
}
  \item{approx_method_SO4}{
    Sets type of interpolation when sulfate is variable (\code{SO4} is a data frame).
    Passed directly to \code{\link{approxfun}}.
    See vignette.
}
  \item{par_key}{
    Special character used to set individual nested elements within some arguments.
    See vignette for examples.
}
  \item{value}{
    Sets the type of output.
    Default (\code{ts}) returns a data frame with the complete time series.
    To exclude those rows before the end of \code{startup}, use \code{tsel}.
    For a summary (after \code{startup} only), use (\code{summ}, \code{sum}, or \code{summary}).
}
  \item{warn}{
    Set to \code{FALSE} to suppress warnings.
}
}
\details{
  The function simulates anaerobic biodegradation.
}
\value{
  By default, a data frame with most input and all output variables in a time series.
  Each row is separated by the value set with \code{delta_t}.
  Typically the primary variable of interest is methane emission, which is returned as a total (g) and rate (g/d), overall or normalized to COD or VS loading.
  These columns all start with \code{CH4}.
  Microbial biomass (g) is given in columns with names that match those used for the names of the groups.
  For more details, see the Output section in the vignette.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Sasha D. Hafner and Frederik R. Dalby
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Simulation 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model with default parameter values and input variables
out1 <- abm()

# Plot cumulative CH4 emission (g) and emission rate (g/d)
plot(CH4_emis_cum ~ time, data = out1, type = 'l')
plot(CH4_emis_rate ~ time, data = out1, type = 'l')

# Plot microbial biomass
matplot(out1$time, out1[, nn <- c('m1', 'm2', 'm3', 'm4', 'm5')],
        type = 'l', lty = 1, xlab = 'Time (d)', ylab = 'Microbial biomass (g)')
legend('topleft', nn, col = 1:6, lty = 1)

# Simulation 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Anaerobic digester with variable temperature
temp_dat <- data.frame(time = 300 + c(0,  1,  6, 100),
                       temp_C =    c(35, 25, 35,  25))

out2 <- abm(600, 1, add_pars = list(temp_C = temp_dat, resid_frac = 0.95,
                                    slurry_mass = 0.95 * 33333, slurry_prod_rate = 500),
            approx_method_temp = 'constant')

# Plot results
plot(temp_C ~ time, type = 'l', col = 'red', data = out2)
abline(v = temp_dat$time, lty = 2, col = 'gray45')

matplot(out2$time, out2[, nn <- c('m1_conc', 'm2_conc', 'm3_conc', 'm4_conc', 'm5_conc')],
        type = 'l', lty = 1, xlab = 'Time (d)', ylab = 'Biomass conc. (g/kg)')
legend('topleft', nn, col = 1:5, lty = 1)
abline(v = temp_dat$time, lty = 2, col = 'gray45')

plot(CH4_emis_rate_VS ~ time, data = out2, type = 'l')
abline(v = temp_dat$time, lty = 2, col = 'gray45')

plot(VFA_conc ~ time, type = 'l', col = 'purple', data = out2)
abline(v = temp_dat$time, lty = 2, col = 'gray45')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
