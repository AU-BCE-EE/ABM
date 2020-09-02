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
mmm(days = 365,                              
    delta_t = 1,                             
    mng_pars = list(slurry_prod_rate = 1000, 
                   slurry_mass = 0,          
                   max_slurry_mass = 33333,  
                   resid_frac = 0.10,        
                   area = 11,                
                   temp_C = 23),
    man_pars = list(conc_fresh = list(S2 = 0.0, SO4 = 0.2, TAN = 1.0, VFA = 4.2, Sp = 65, COD = 160),
                   pH = 7), # Note list not c() so SO4 can be data frame
    grp_pars = list(yield = c(default = 0.04, sr1 = 0.065),
                    xa_fresh = c(default = 0.001, sr1 = 0.001),
                    xa_init = c(m1 = 0.01, m2 = 0.01, m3 = 0.01, p1 = 0.01, p2 = 0.01, sr1 = 0.01),
                    decay_rate = c(m1 = 0.02, m2 = 0.02, m3 = 0.02, p1 = 0.02, p2 = 0.02, sr1 = 0.02),
                    ks_coefficient = c(m1 = 0.5, m2 = 1.5, m3 = 1.0, p1 = 1.0, p2 = 1.0, sr1 = 0.4),
                    #resid_enrich = c(m1 = 0.5, m2 = 0.5, m3 = 0.5, p1 = 0.5, p2 = 0.5, sr1 = 0.5),
                    resid_enrich = c(m1 = 0.0, m2 = 0.0, m3 = 0.0, p1 = 0.0, p2 = 0.0, sr1 = 0.0),
                    qhat_opt = c(m1 = 8, m2 = 13.33, m3 = 5.75, p1 = 2.77, p2 = 0.72, sr1 = 8.3),    
                    T_opt = c(m1 = 313, m2 = 313, m3 = 303, p1 = 293, p2 = 283, sr1 = 313),
                    T_min = c(m1 = 295.31, m2 = 295.31, m3 = 285.31, p1 = 275.31, p2 = 265.31, sr1 = 273),
                    T_max = c(m1 = 320.67, m2 = 320.67, m3 = 310.67, p1 = 300.67, p2 = 290.67, sr1 = 320.67),
                    ki_NH3_min = c(m1 = 0.01, m2 = 0.015, m3 = 0.015, p1 = 0.015, p2 = 0.015, sr1 = 0.015),
                    ki_NH3_max = c(m1 = 0.10, m2 = 0.131, m3 = 0.131, p1 = 0.131, p2 = 0.131, sr1 = 0.131),
                    ki_NH4_min = c(m1 = 1.70, m2 = 2.714, m3 = 2.714, p1 = 2.714, p2 = 2.714, sr1 = 2.714),
                    ki_NH4_max = c(m1 = 3.10, m2 = 4.764, m3 = 4.764, p1 = 4.764, p2 = 4.764, sr1 = 4.764),
                    pH_upr = c(m1 = 8.0, m2 = 8.0, m3 = 8.0, p1 = 8.0, p2 = 8.0, sr1 = 8.0),
                    pH_lwr = c(m1 = 6.5, m2 = 6.0, m3 = 6.5, p1 = 6.5, p2 = 6.5, sr1 = 6.0)),
    mic_pars = list(ks_SO4 = 0.0067,
                    ki_H2S_meth = 0.23,
                    ki_H2S_sr = 0.25,
                    alpha_opt = 0.015,
                    alpha_T_opt = 313,
                    alpha_T_min = 273,
                    alpha_T_max = 320.67),
    chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.57, CO2_aer = 1.3, CO2_sr = 1.3), kl = c(H2S = 0.032, oxygen = 0.415)),  
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
    List of management parameters.
    % Should list elements
    % Add details here and for other parameter arguments
    Each element may be a length-one vector or, in some cases, a list.
    See usage, examples, and vignettes.
}
  \item{man_pars}{
    List of manure parameters.
}
  \item{grp_pars}{
    List containing parameter values for the microbial groups included in the model.
}
  \item{mic_pars}{
    List of other microbial parameters.
}
  \item{chem_pars}{
    List of other microbial parameters.
}
  \item{add_pars}{
    List of other additional parameters not set in other arguments, or to be modified from defaults.
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
%%     ~~Describe \code{approx_method_temp} here~~
}
  \item{approx_method_pH}{
%%     ~~Describe \code{approx_method_pH} here~~
}
  \item{approx_method_SO4}{
%%     ~~Describe \code{approx_method_SO4} here~~
}
  \item{par_key}{
%%     ~~Describe \code{par_key} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{warn}{
%%     ~~Describe \code{warn} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
