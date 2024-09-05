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
    times = NULL,
    wthr_pars = ABM::wthr_pars2.0,
    evap_pars = list(evap = 0.5 * et(temp_C = ABM::wthr_pars2.0$temp_air_C, 
                                     pres_kpa = ABM::wthr_pars2.0$pres_kpa, 
                                     rs = ABM::wthr_pars2.0$rs)), 
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
                    cover = 'none',
                    resid_enrich = 0.9,
                    slopes = c(urea = NA, slurry_prod_rate = NA),
                    graze = c(start = 'May', duration = 0, hours_day = 0),
                    scale = c(ks_coefficient = 1, qhat_opt = 1, 
                              xa_fresh = 1, yield = 1, alpha_opt = 1)),
    man_pars = ABM::man_pars2.0,
    init_pars = list(conc_init = man_pars$conc_fresh),
    grp_pars = ABM::grp_pars2.0,
    mic_pars = ABM::mic_pars2.0,
    chem_pars = ABM::chem_pars2.0,
    arrh_pars = ABM::arrh_pars2.0,
    resp = TRUE,
    add_pars = NULL,
    pars = NULL,
    startup = 0,
    starting = NULL,
    approx_method_temp = 'linear', 
    approx_method_pH = 'linear',
    approx_method_vent_air = 'linear',
    par_key = '\\\\.',
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
  \item{times}{
    Alternative to \code{delta_t}, use to return specific times only.
    Only available with variable slurry mass.
}
  \item{wthr_pars}{
    Weather inputs.
    List with elements `temp_air_C` in degree C, `RH` relative humidity in percent, `rain` is precipitation in kg/m2/day, `pres_kpa` is ambient pressure in kPa, `rs` is solar radiation in MJ/m2/d.
    See \code{\link{parameters}}.
}
  \item{evap_pars}{
    Evaporation inputs.
}
  \item{mng_pars}{
    List of management parameters including slurry production rate (kg/day) (\code{slurry_prod_rate}),
    initial slurry mass in the system (kg) (\code{slurry_mass}). If (\code{slurry_mass}) is a data.frame with column names "time", "slurry_mass", then (\code{slurry_prod_rate}), (\code{resid_depth}), (\code{empty_int}), (\code{wash_water}), (\code{wash_int}) are calculated from the (\code{slurry_mass}) data.frame.
    depth of slurry storage (m) (\code{storage_depth}),
    height of slurry in the system after being emptied (m) (\code{resid_depth}),
    surface area of the barn floor that is accesible by the animals (m^2) (\code{floor_area}),
    surface area of the slurry storage (m^2) (\code{area}),
    removal frequency of the slurry in the slurry storage (days), where a value of 0 = "never remove slurry" (\code{empty_int}),
    slurry temperature (deg. C) (\code{temp_C}). The model provides some datasets for outside slurry temperatures in Denmark \code{\link{temperatures}},
    wash water used per washing event (kg) (\code{wash_water}),
    time between washing events (days), where NA = "washing not applied" (\code{wash_int}),
    time that the animal section is empty after a batch of animals has been finished (days) (\code{rest_d}), argument for setting the type of cover on an outside   storage tank (\code{cover}). The type of cover only affects the ammonia emission factor, with emission factors and types of covers taken from Sommer et al. 2022, e.g. cover = \code{none} means no cover, cover = \code{tent} reduces ammonia emission by 80 percent, cover = \code{straw} reduces ammonia emission by 33 percent.
    Factor that determines the degree of enrichment of microbes and organic matter components in the residual slurry after emptying (\code{resid_enrich}),
    and \code{slopes} is an expert user argument for changing the excretion rate of urea or slurry production during a growth period of an animal. 
    By default this argument is not used (\code{slopes}),
    (\code{graze}) is a list with parameters to adjust the grazing activity outdoor, which affects the amount of slurry excreted in the barn. 
    \code{duration} is the number of days per year that the animals are allowed to graze and \code{hours_day} are the number of hours animals graze per day during the grazing period, the \code{start} argument is the month, where grazing starts, e.g. "May". When the (\code{graze}) arguments are not 0, then the model simulation is with time reference to 1. January 00:00:00. Effectively, if start = \code{May}, the model simulates grazing from time = 120. 
    Therefore, when combined with variable temperature, temperature should also be with time = 0 referencing to 1 January 00:00:00. 
    \code{scale} is a multiplier for adjusting some \code{grp_pars} parmeters for all microbes at the same time (\code{scale}).
    Each element may be a length-one vector or, in some cases, a list.
    See usage, examples, and vignette.
}
  \item{man_pars}{
    List of manure parameters. The default and alternative parameter values are detailed in \code{\link{parameters}}.
    Elements are:
    \code{conc_fresh}, the concentration of slurry components in fresh slurry
    \code{pH}, pH of the slurry entering the system and in the slurry in the system
    The \code{conc_fresh} element is a list itself, with the following elements:
    \code{sulfide} (sulfide), \code{urea} (urea), \code{sulfate} (sulfate), \code{TAN} (ammonia nitrogen), \code{starch} (starch and sugars),
    \code{VFA} (volatile fatty acids), \code{xa_dead} (dead microbial populations), \code{Cfat} (crude fat),
    \code{CP} (crude protein), \code{RFd} (residual fiber that is degradable: defined as neutral detergent fiber - indigestible neutral detergent fiber), \code{iNDF} (indigestible neutral detergent fiber), \code{VSd} (degradable, volatile solids), \code{VSd_A} (degradable volatile solids for Arrhenius calculation), \code{VSnd_A} (non-degradable volatile solids for Arrhenius calcualtion), \code{ash} (ash content)
    Concentrations are on an elemental basis (N/kg wet manure and S/kg wet manure) or else gCOD/kg wet manure (organic components). However, for VSd_A, VSnd_A and ash the unit is g/kg wet manure.
}
  \item{init_pars}{
    List of parameters that define the intial concentrations of organic matter components in the slurry.
    The elements are similar to the \code{conc_fresh} elements in 'man_pars'. 
    By default the element values are the same as values in \code{conc_fresh}
}

  \item{grp_pars}{
   List containing parameter values for the microbial groups included in the model. m designates methanogen and sr designates sulfate reducer. The default and alternative parameter values are detailed in \code{\link{parameters}}.
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
    List of other microbial parameters. The default and alternative parameter values are detailed in \code{\link{parameters}}.
    Elements are:
    \code{ks_SO4}, Half sulfate saturation constant for sulfate reducers (g sulfate-sulfur/kg);
    \code{km_urea}, half urea-N saturation constant for Michaelis-Menten urea hydrolysis (g urea-N/kg).
}
  \item{chem_pars}{
    List of other chemical parameters. 
    The default parameter values are detailed in \code{\link{parameters}}.
    Elements are:
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
    \code{CO2_ureo}, carbon dioxide productivity coefficient for urea hydrolysis (g CO2/g urea-N);
    \code{N_CP}, nitrogen per COD of crude protein (g N/g CP-COD);
    \code{C_N_urea}, carbon to nitrogen (g C/g urea-N);
    Other conversion factors designated \code{C_} are carbon to COD of the named substrate (gC/gCOD) 
}

  \item{arrh_pars}{
    List of hydrolysis parameters that are related to Arrhenius kinetics. 
    The default parameter values are detailed in \code{\link{parameters}}.
    Elements are:
    \code{lnA}, frequency factor used in danish methane emission inventory calculations, only used for comparing to ABM model (g CH4/ (kgVS * hour));
    \code{E_CH4}, Activation energy used in danish methane emission inventory calculations, only used for comparing to ABM model (J/mol);
    \code{R}, universal gas constant (J/(K * mol));
    \code{VS_CH4}, volatile solids consumed per methane produced (g VS/g CH4); 
    \code{A}, Arrhenius factor for calculating the temperature sensitive hydrolysis rate constants of organic matter components (1/day);
    \code{E}, Activation energy for calculating the temperature sensitive hydrolysis rate constant of organic matter components (J/mol);
    \code{kl}, mass transfer coefficient of ammonia from the slurry pit (m/day) and from the barn floor (m/day) and mass transfer coefficient of H2S from the slurry (m/day).

}
 \item{resp}{
    Set to \code{TRUE} (default) to include surface respiration. 
    Can be \code{FALSE} for e.g. anaerobic digestion simulation.
}
    
 \item{pH_inhib_overrule}{
    Default is \code{FALSE}. 
    Can be \code{TRUE} if the inhibition should be governed by a pH inhibition function, rather than the default where inhibition is governed by concentrations of NH3, NH4, H2S and protonated VFA.
}
        
  \item{add_pars}{
    List of other additional parameters not set in other arguments, or to be modified from defaults.
    See vignette.
}
  \item{pars}{
    Combined parameters, as an alternative to all the `*_pars` arguments.
    Used internally when \code{startup} > 0 and may be convenient to use in an \code{abm} call in some cases.
}
  \item{startup}{
    Number of startup times the simulation should be run before running a final time.
    Results from one simulation become starting conditions for next.
    Used to spin up or stabilize microbial population and substrates so that initial conditions do not affect results.
    Length-one integer vector.
    Optional.
}
  \item{starting}{
    Starting conditions.
    Output (\code{value = 'ts'}) from a previous \code{abm()} call.
    Values here for microbial biomass and slurry mass override those set in other arguments.
    Used to start one simulation where another ended.
    Optional.
}
  \item{approx_method_temp}{
    Sets type of interpolation when temperature is variable (\code{temp_c} is a data frame).
    Passed to \code{\link{approxfun}} through hidden (not exported) \code{makeTimeFunc()}.
    See vignette. Suggested temperature sets for outside slurry storage are provided in \code{\link{temperatures}}.
}
  \item{approx_method_pH}{
    Sets type of interpolation when pH is variable (\code{pH} is a data frame).
    Passed to \code{\link{approxfun}} through hidden (not exported) \code{makeTimeFunc()}.
    See vignette.
}
  \item{approx_method_vent_air}{
    Sets type of interpolation when ventilation flow rate is variable.
    Passed to \code{\link{approxfun}} through hidden (not exported) \code{makeTimeFunc()}.
    See vignette.
}
  \item{par_key}{
    Special character used to set individual nested elements within some arguments.
    See vignette for examples.
}
  \item{value}{
    Sets the type of output.
    Default (\code{'ts'}) returns a data frame with the complete time series.
    For a summary, use \code{'summ'} (or \code{'sum'} or \code{'summary'}).
    code{'eff'} will return only those rows where slurry was removed--this option is convenient for working with effluent slurry.
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
Dalby, F.R., Hafner, S.D., Petersen, S.O., Vanderzaag, A., Habtewold, J., Dunfield, K., Chantigny, M.H., Sommer, S.G., 2021. A mechanistic model of methane emission from animal slurry with a focus on microbial groups. PLOS ONE 16, e0252881. \doi{10.1371/journal.pone.0252881}

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. \doi{10.1021/acs.est.2c08891}

Dalby, F.R., Ambrose, H.W.,  Poulsen, J.S., Nielsen, J.L., Adamsen, A.P.S., 2023. Pig slurry organic matter transformation and methanogenesis at ambient storage temperatures. JEQ. \doi{10.1002/jeq2.20512}

Dalby, F.R., Hafner, S.D., Ambrose, H.W., Adamsen, A.P.S., 2024. Pig manure degradation and carbon emission: Measuring and modeling combined aerobic–anaerobic transformations. JEQ. \doi{10.1002/jeq2.20603}

Sommer, S.G, Hafner S.D., Laubach, J., van der Weerden, T.J., Pacholski, A., 2022. Model for calculating ammonia emission from stored animal liquid manure. 
Biosyst Eng \doi{10.1016/j.biosystemseng.2022.08.007}
} 
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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
