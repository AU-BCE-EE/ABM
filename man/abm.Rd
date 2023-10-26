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
    mng_pars = list(slurry_prod_rate = 5700,   
                  slurry_mass = 39000,          
                  storage_depth = 0.6,        
                  resid_depth = 0.05,         
                  floor_area = 650,           
                  area = 715,                 
                  empty_int = 42,          
                  temp_C = 20,
                  wash_water = 75000,            
                  wash_int = NA,
                  rest_d = 5,
                  cover = NA,
                  resid_enrich = 0.9,
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  graze = c(start = 'May', duration = 0, hours_day = 0),
                  scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1)),
    grz_pars = grz_pars2.0,
    man_pars = man_pars2.0,
    init_pars = list(conc_init = man_pars$conc_fresh),
    grp_pars = grp_pars2.0,
    mic_pars = mic_pars2.0,
    chem_pars = chem_pars2.0,
    arrh_pars = arrh_pars2.0,
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
    depth of slurry storage (m) (\code{storage_depth}),
    height of slurry in the system after being emptied (m) (\code{resid_depth}),
    surface area of the barn floor that is accesible by the animals (m^2) (\code{floor_area}),
    surface area of the slurry storage (m^2) (\code{area}),
    removal frequency of the slurry in the slurry storage (days), where a value of 0 = "never remove slurry"" (\code{empty_int}),
    slurry temperature (deg. C) (\code{temp_C}),
    wash water used per washing event (kg) (\code{wash_water}),
    time between washing events (days), where NA = "washing not applied" (\code{wash_int}),
    time that the animal section is empty after a batch of animals has been finished (days) (\code{rest_d}),
    factor that determines the degree of enrichment of microbes and organic matter components in the residual slurry after emptying   (\code{resid_enrich}),
    slopes is an expert user argument for changing the increasing the excretion rate of urea or slurry production during a growth period of an animal     . By default this argument is not used (\code{slopes}),
    multiplyer arguments for adjusting some grp_pars arguments for all microbes at the same time (\code{scale}),
    
    Each element may be a length-one vector or, in some cases, a list.
    See usage, examples, and vignettes.
}
  \item{man_pars}{
    List of manure parameters.
    Elements are:
    \code{conc_fresh}, the concentration of slurry components in fresh slurry;
    \code{pH}, pH of the slurry entering the system and in the slurry in the system.
    The \code{conc_fresh} element is a list itself, with the following elements:
    \code{sulfide} (sulfide), \code{urea} (urea), \code{sulfate} (sulfate), \code{TAN} (ammonia nitrogen), \code{starch} (starch and sugars),
    \code{VFA} (volatile fatty acids), \code{xa_dead} (dead microbial populations), \code{Cfat} (crude fat),
    \code{CP} (crude protein), \code{RFd} (residual fiber that is degradable: defined as neutral detergent fiber - indigestible neutral detergent fiber), \code{iNDF} (indigestible neutral detergent fiber), \code{VSd} (degradable, volatile solids), \code{VSd_A} (degradable volatile solids for Arrhenius calculation), \code{VSnd_A} (non-degradable volatile solids for Arrhenius calcualtion), \code{ash} (ash content)
    Concentrations are on an elemental basis (N/kg wet manure and S/kg wet manure) or else gCOD/kg wet manure (organic components). However, for VSd_A, VSnd_A and ash the unit is g/kg wet manure.
}
  \item{init_pars}{
    List of parameters that define the intial concentrations of organic matter components in the slurry.
    The elements are similar to the \code{conc_fresh} elements in 'man_pars'. By default the element values are the same as values in \code{conc_fresh}
    \code{graze_start}, month in which grazing (for cattle) starts, e.g. in "may"" (month);
    \code{graze_days}, days per year that the animals are outside grazing (days/year);
    \code{graze_hours}, hours per day the animals are outside grazing (hours/day).
}

  \item{grp_pars}{
   List containing parameter values for the microbial groups included in the model. m designates methanogen and sr designates sulfate reducer.
   Elements are:
   \code{yield}, Growth yield of microbial groups (g COD-biomass / g COD-substrate) ;
   \code{xa_fresh}, Concentration of microbial groups in the slurry entering the system (g COD-biomass/kg);
   \code{xa_init}, Initial concentration of microbial groups in the slurry (g COD-biomass/kg);
   \code{decay_rate}, Deacay rate of microbial groups (1/day);
   \code{ks_coefficient}, Coefficient used to tune the half saturation constant (g COD-substrate/kg);
   \code{qhat_opt}, Maximum specific substrate utilization rate of microbial groups at their optimum growth temperature (g COD-substrate/g COD-biomass/day);
   \code{T_opt}, Optimum growth temperature of microbial groups (K);
   \code{T_min}, Minimum viable growth temperature of microbial groups (K);
   \code{T_max}, Maximum vaible growth temperature of microbial groups (K);
   \code{ki_NH3_min}, NH3 concentration at which NH3 inhibition starts (g TAN-nitrogen/kg);
   \code{ki_NH3_max}, NH3 concentration wher microbial groups are 100 percent inhibited (g TAN-nitrogen/kg);
   \code{ki_NH4_min}, NH4 concentration at which NH3 inhibition starts (g TAN-nitrogen/kg);
   \code{ki_NH4_max}, NH4 concentration wher microbial groups are 100 percent inhibited (g TAN-nitrogen/kg);
   \code{ki_H2S_slope}, inhibition parameter1 for calculating H2S inhibition factor (g H2S-sulfur/(kg pr pH unit));
   \code{ki_H2S_int}, inhibition parameter2 for calculating H2S inhibition factor (g H2S-sulfur/kg);
   \code{ki_H2S_min}, inhibition parameter3 for calculating H2S inhibition factor (g H2S-sulfur/kg);

   Each of these elements is a vector with a named element for each group.
   See usage, examples, and vignettes.
}
  \item{mic_pars}{
    List of other microbial parameters.
    Elements are:
    \code{ks_SO4}, Half sulfate saturation constant for sulfate reducers (g sulfate-sulfur/kg);
    \code{km_urea}, half urea-N saturation constant for Michaelis-Menten urea hydrolysis (g urea-N/kg);
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
    \code{CH4}, Methane productivity coefficient (g COD-substrate/g CH4) ;
    \code{xa_dead}, COD to g dry matter mass of dead microbial populations (g COD/g dry xa dead);
    \code{RFd}, COD to g RFd (g COD/g RFd);
    \code{iNDF}, COD to g iNDF (g COD/g iNDF);
    \code{starch}, COD to g starch (g COD/g starch);
    \code{Cfat}, COD to g Cfat (g COD/g Cfat);
    \code{CP}, COD to g CP (g COD/g CP);
    \code{VFA}, COD to g VFA (g COD/g VFA);
    \code{RFd}, COD to g RFd (g COD/g RFd);
    \code{S}, COD-substrate to sulfur (g sulfur/g COD-substrate) ;
    \code{VS}, COD to g volatile solids (g COD/g volatile solids);
    \code{CO2_aer}, Aerobic carbon dioxide productivity coefficient (g CO2/g COD-substrate);
    \code{CO2_sr}, carbon dioxide productivity coefficient for sulfate reduction (g CO2/g COD-substrate);
    \code{CO2_ureo}, carbon dioxide productivity coefficient for urea hydrolysis (g CO2/g urea-N);
    \code{N_CP}, nitrogen per COD of crude protein (g N/g CP-COD).
    \code{C_N_urea}, carbon to nitrogen (g C/g urea-N)
    \code{kl}, mass transfer coefficient of ammonia from the slurry pit [NH3] (m/day) and from the barn floor [NH3_floor] (m/day) and mass transfer coefficient of H2S from the slurry (m/day).
    Other conversion factors designated \code{C_} are carbon to COD of the named substrate (gC/gCOD) 
}

  \item{arrh_pars}{
    List of hydrolysis parameters that are related to arrhenius kinetics.
    Elements are:
    \code{lnA}, frequency factor used in danish methane emission inventory calculations, only used for comparing to ABM model (g CH4/ (kgVS * hour));
    \code{E_CH4}, Activation energy used in danish methane emission inventory calculations, only used for comparing to ABM model (J/mol);
    \code{R}, universal gas constant (J/(K * mol));
    \code{VS_CH4}, volatile solids consumed per methane produced (g VS/g CH4); 
    \code{A}, Arrhenius factor for calculating the temperature sensitive hydrolysis rate constant of organic amtter compoenents (1/day);
    \code{E}, Activation energy for calculating the temperature sensitive hydrolysis rate constant of organic amtter compoenents (J/mol).

}

  \item{grz_pars}{
    List of parameters used for simualting periods of grazing.
    Elements are:
    \code{graze_start}, month in which grazing (for cattle) starts, e.g. in "may"" (month);
    \code{graze_days}, days per year that the animals are outside grazing (days/year);
    \code{graze_hours}, hours per day the animals are outside grazing (hours/day).
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
matplot(out1$time, out1[, nn <- c('m0', 'm1', 'm2', 'sr1')],
        type = 'l', lty = 1, xlab = 'Time (d)', ylab = 'Microbial biomass (g)')
legend('topleft', nn, col = 1:4, lty = 1)

# Simulation 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model with variable temperature and constant slurry mass
temp_dat <- data.frame(time = 300 + c(0,  1,  6, 100),
                       temp_C =    c(10, 25, 10, 25))

out2 <- abm(600, 1, add_pars = list(temp_C = temp_dat, slurry_mass = 8000, slurry_prod_rate = 0,
            approx_method_temp = 'linear')

# Plot results
plot(temp_C ~ time, type = 'l', col = 'red', data = out2)
abline(v = temp_dat$time, lty = 2, col = 'gray45')

matplot(out2$time, out2[, nn <- c('m0_conc', 'm1_conc', 'm2_conc')],
        type = 'l', lty = 1, xlab = 'Time (d)', ylab = 'Biomass conc. (g/kg)')
legend('topleft', nn, col = 1:3, lty = 1)
abline(v = temp_dat$time, lty = 2, col = 'gray45')

plot(CH4_emis_rate_VS ~ time, data = out2, type = 'l')
abline(v = temp_dat$time, lty = 2, col = 'gray45')

plot(VFA_conc ~ time, type = 'l', col = 'purple', data = out2)
abline(v = temp_dat$time, lty = 2, col = 'gray45')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
