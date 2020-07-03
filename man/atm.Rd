\name{mmm}
\alias{mmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mmm(days = 365, delta_t = 1, mng_pars = list(slurry_prod_rate = 1000, slurry_mass = 0, max_slurry_mass = 33333, resid_frac = 0.95, area = 11, temp_C = 23), man_pars = list(conc_fresh = list(S2 = 0, SO4 = 0.2, TAN = 1, VFA = 4.2, Sp = 63.7, COD = 163.6), pH = 7), grp_pars = list(yield = c(m1 = 0.04, m2 = 0.04, m3 = 0.04, p1 = 0.04, p2 = 0.04, sr1 = 0.065), xa_fresh = c(m1 = 0.001, m2 = 0.01, m3 = 0.001, p1 = 0.001, p2 = 0.001, sr1 = 0.001), xa_init = c(m1 = 0.01, m2 = 0.01, m3 = 0.01, p1 = 0.01, p2 = 0.01, sr1 = 0.01), decay_rate = c(m1 = 0.02, m2 = 0.02, m3 = 0.02, p1 = 0.02, p2 = 0.02, sr1 = 0.02), ks_coefficient = c(m1 = 0.5, m2 = 1.5, m3 = 1, p1 = 1, p2 = 1, sr1 = 0.4), resid_enrich = c(m1 = 2, m2 = 2, m3 = 2, p1 = 2, p2 = 2, sr1 = 2), qhat_opt = c(m1 = 8, m2 = 13.33, 
    m3 = 5.75, p1 = 2.77, p2 = 0.72, sr1 = 8.3), T_opt = c(m1 = 313, m2 = 313, m3 = 303, p1 = 293, p2 = 283, sr1 = 313), T_min = c(m1 = 295.31, m2 = 295.31, m3 = 285.31, p1 = 275.31, p2 = 265.31, sr1 = 273), T_max = c(m1 = 320.67, m2 = 320.67, m3 = 310.67, p1 = 300.67, p2 = 290.67, sr1 = 320.67), ki_NH3_min = c(m1 = 0.01, m2 = 0.015, m3 = 0.015, p1 = 0.015, p2 = 0.015, sr1 = 0.015), ki_NH3_max = c(m1 = 0.1, m2 = 0.131, m3 = 0.131, p1 = 0.131, p2 = 0.131, sr1 = 0.131), ki_NH4_min = c(m1 = 1.7, m2 = 2.714, 
    m3 = 2.714, p1 = 2.714, p2 = 2.714, sr1 = 2.714), ki_NH4_max = c(m1 = 3.1, m2 = 4.764, m3 = 4.764, p1 = 4.764, p2 = 4.764, sr1 = 4.764), pH_upr = c(m1 = 8, m2 = 8, m3 = 8, p1 = 8, p2 = 8, sr1 = 8), pH_lwr = c(m1 = 6.5, m2 = 6, m3 = 6.5, p1 = 6.5, p2 = 6.5, sr1 = 6)), mic_pars = list(ks_SO4 = 0.0067, ki_H2S_meth = 0.23, ki_H2S_sr = 0.25, alpha_opt = 0.015, alpha_T_opt = 313, alpha_T_min = 273, alpha_T_max = 320.67), chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.57, CO2_aer = 1.3, CO2_sr = 1.3), kl = c(H2S = 0.032, oxygen = 0.415)), add_pars = NULL, startup = -Inf, starting = NULL, approx_method_temp = "linear", approx_method_pH = "linear", approx_method_SO4 = "linear", par_key = "\\.", value = "ts", warn = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{days}{
%%     ~~Describe \code{days} here~~
}
  \item{delta_t}{
%%     ~~Describe \code{delta_t} here~~
}
  \item{mng_pars}{
%%     ~~Describe \code{mng_pars} here~~
}
  \item{man_pars}{
%%     ~~Describe \code{man_pars} here~~
}
  \item{grp_pars}{
%%     ~~Describe \code{grp_pars} here~~
}
  \item{mic_pars}{
%%     ~~Describe \code{mic_pars} here~~
}
  \item{chem_pars}{
%%     ~~Describe \code{chem_pars} here~~
}
  \item{add_pars}{
%%     ~~Describe \code{add_pars} here~~
}
  \item{startup}{
%%     ~~Describe \code{startup} here~~
}
  \item{starting}{
%%     ~~Describe \code{starting} here~~
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
%%  ~~who you are~~
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

## The function is currently defined as
function (days = 365, delta_t = 1, mng_pars = list(slurry_prod_rate = 1000, 
    slurry_mass = 0, max_slurry_mass = 33333, resid_frac = 0.95, 
    area = 11, temp_C = 23), man_pars = list(conc_fresh = list(S2 = 0, 
    SO4 = 0.2, TAN = 1, VFA = 4.2, Sp = 63.7, COD = 163.6), pH = 7), 
    grp_pars = list(yield = c(m1 = 0.04, m2 = 0.04, m3 = 0.04, 
        p1 = 0.04, p2 = 0.04, sr1 = 0.065), xa_fresh = c(m1 = 0.001, 
        m2 = 0.01, m3 = 0.001, p1 = 0.001, p2 = 0.001, sr1 = 0.001), 
        xa_init = c(m1 = 0.01, m2 = 0.01, m3 = 0.01, p1 = 0.01, 
            p2 = 0.01, sr1 = 0.01), decay_rate = c(m1 = 0.02, 
            m2 = 0.02, m3 = 0.02, p1 = 0.02, p2 = 0.02, sr1 = 0.02), 
        ks_coefficient = c(m1 = 0.5, m2 = 1.5, m3 = 1, p1 = 1, 
            p2 = 1, sr1 = 0.4), resid_enrich = c(m1 = 2, m2 = 2, 
            m3 = 2, p1 = 2, p2 = 2, sr1 = 2), qhat_opt = c(m1 = 8, 
            m2 = 13.33, m3 = 5.75, p1 = 2.77, p2 = 0.72, sr1 = 8.3), 
        T_opt = c(m1 = 313, m2 = 313, m3 = 303, p1 = 293, p2 = 283, 
            sr1 = 313), T_min = c(m1 = 295.31, m2 = 295.31, m3 = 285.31, 
            p1 = 275.31, p2 = 265.31, sr1 = 273), T_max = c(m1 = 320.67, 
            m2 = 320.67, m3 = 310.67, p1 = 300.67, p2 = 290.67, 
            sr1 = 320.67), ki_NH3_min = c(m1 = 0.01, m2 = 0.015, 
            m3 = 0.015, p1 = 0.015, p2 = 0.015, sr1 = 0.015), 
        ki_NH3_max = c(m1 = 0.1, m2 = 0.131, m3 = 0.131, p1 = 0.131, 
            p2 = 0.131, sr1 = 0.131), ki_NH4_min = c(m1 = 1.7, 
            m2 = 2.714, m3 = 2.714, p1 = 2.714, p2 = 2.714, sr1 = 2.714), 
        ki_NH4_max = c(m1 = 3.1, m2 = 4.764, m3 = 4.764, p1 = 4.764, 
            p2 = 4.764, sr1 = 4.764), pH_upr = c(m1 = 8, m2 = 8, 
            m3 = 8, p1 = 8, p2 = 8, sr1 = 8), pH_lwr = c(m1 = 6.5, 
            m2 = 6, m3 = 6.5, p1 = 6.5, p2 = 6.5, sr1 = 6)), 
    mic_pars = list(ks_SO4 = 0.0067, ki_H2S_meth = 0.23, ki_H2S_sr = 0.25, 
        alpha_opt = 0.015, alpha_T_opt = 313, alpha_T_min = 273, 
        alpha_T_max = 320.67), chem_pars = list(COD_conv = c(CH4 = 0.2507, 
        S = 0.5015, VS = 0.69, CO2_anaer = 0.57, CO2_aer = 1.3, 
        CO2_sr = 1.3), kl = c(H2S = 0.032, oxygen = 0.415)), 
    add_pars = NULL, startup = -Inf, starting = NULL, approx_method_temp = "linear", 
    approx_method_pH = "linear", approx_method_SO4 = "linear", 
    par_key = "\\.", value = "ts", warn = TRUE) 
{
    if (!is.null(starting) & is.data.frame(starting)) {
        message("Using starting conditions from `starting` argument")
        grp_pars[["xa_init"]] <- starting[nrow(starting), paste0(names(grp_pars[["qhat_opt"]]), 
            "_conc")]
        names(grp_pars[["xa_init"]]) <- names(grp_pars[["qhat_opt"]])
    }
    pars <- c(mng_pars, man_pars, grp_pars, mic_pars, chem_pars)
    if (!is.null(add_pars) && length(add_pars) > 0 && any(ii <- grepl(par_key, 
        names(add_pars)))) {
        pe.pars <- add_pars[ii]
        sa.pars <- add_pars[!ii]
        apnames <- names(pe.pars)
        pe.pars <- as.numeric(pe.pars)
        split.pars <- strsplit(apnames, par_key)
        pnames <- sapply(split.pars, "[[", 1)
        enames <- sapply(split.pars, "[[", 2)
        names(pe.pars) <- enames
        pe.pars <- split(pe.pars, pnames)
        add_pars <- c(sa.pars, pe.pars)
    }
    if (!is.null(add_pars) && length(add_pars) > 0) {
        if (any(bad.names <- !names(add_pars) \%in\% names(pars))) {
            stop("Some `add_pars` names not recognized as valid parameters: ", 
                names(add_pars)[bad.names])
        }
        for (i in names(add_pars)) {
            if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 
                1) {
                pars[[i]][names(add_pars[[i]])] <- add_pars[[i]]
            }
            else {
                pars[[i]] <- add_pars[[i]]
            }
        }
    }
    grp_nms <- names(pars$qhat_opt)
    for (i in names(grp_pars)) {
        ppo <- pars[[i]]
        p_nms <- names(pars[[i]])
        if (any(p_nms == "default")) {
            pars[[i]][grp_nms] <- pars[[i]]["default"]
            if (any(p_nms != "default")) {
                pars[[i]][p_nms[p_nms != "default"]] <- ppo[p_nms[p_nms != 
                  "default"]]
            }
        }
        if (any(p_nms == "all")) {
            pars[[i]][grp_nms] <- pars[[i]]["all"]
        }
        pars[[i]] <- pars[[i]][grp_nms]
    }
    if (!all.equal(names(pars$yield), names(pars$xa_fresh), names(pars$xa_init), 
        names(pars$decay_rate), names(pars$ks_coefficient), names(pars$resid_enrich), 
        names(pars$qhat_opt), names(pars$T_opt), names(pars$T_min), 
        names(pars$T_max), names(pars$ki_NH3_min), names(pars$ki_NH3_max), 
        names(pars$ki_NH4_min), names(pars$ki_NH4_max), names(pars$pH_lwr), 
        names(pars$upr))) {
        stop("Microbial groups, i.e., element names in `grp_pars`, must match.")
    }
    if (is.data.frame(pars$temp_C)) {
        temp <- pars$temp_C$temp_C
        ttime <- pars$temp_C$time
        temp_C_fun <<- approxfun(ttime, temp, method = approx_method_temp, 
            yleft = temp[1], yright = temp[length(temp)], rule = 2)
    }
    else {
        temp_C_fun <<- function(x) return(pars$temp_C)
    }
    if (is.data.frame(pars$pH)) {
        tpH <- pars$pH$pH
        ttime <- pars$pH$time
        pH_fun <<- approxfun(ttime, tpH, method = approx_method_pH, 
            yleft = tpH[1], yright = tpH[length(tpH)], rule = 2)
    }
    else {
        pH_fun <<- function(x) return(pars$pH)
    }
    if (is.data.frame(pars$conc_fresh[["SO4"]])) {
        tSO4 <- pars$conc_fres$SO4$SO4
        ttime <- pars$conc_fres$SO4$time
        SO4_fun <<- approxfun(ttime, tSO4, method = approx_method_SO4, 
            yleft = tSO4[1], yright = tSO4[length(tSO4)], rule = 2)
    }
    else {
        SO4_fun <<- function(x) return(pars$conc_fresh[["SO4"]])
    }
    if (is.numeric(pars$slurry_mass)) {
        dat <- mmm_regular(days = days, delta_t = delta_t, pars = pars)
    }
    else if (is.data.frame(pars$slurry_mass)) {
        dat <- mmm_variable(days = days, delta_t = delta_t, pars = pars, 
            warn = warn)
    }
    mic_names <- names(pars$qhat)
    conc.names <- c("NH4", "NH3", "Sp", "VFA", "sulfide", "sulfate", 
        mic_names)
    dat_conc <- dat[, conc.names]/dat$slurry_mass
    names(dat_conc) <- paste0(names(dat_conc), "_conc")
    dat <- cbind(dat, dat_conc)
    dat$temp_C <- temp_C_fun(dat$time)
    if (is.numeric(pars$pH) | is.data.frame(pars$pH)) {
        dat$pH <- pH_fun(dat$time)
    }
    else if (pars$pH == "calc") {
        dat$pH <- H2SO4_titrat(dat$sulfate_conc)
    }
    else {
        stop("Problem with pH input (bee721)")
    }
    conc_fresh <- unlist(pars$conc_fresh)
    names(conc_fresh) <- paste0(names(conc_fresh), "_conc_fresh")
    dat <- cbind(dat, t(conc_fresh))
    dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
    dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
    dat$CH4_emis_rate_slurry <- dat$CH4_emis_rate/(dat$slurry_mass/1000)
    dat$CH4_flux <- dat$CH4_emis_rate/pars$area
    dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
    dat$CO2_emis_rate_slurry <- dat$CO2_emis_rate/(dat$slurry_mass/1000)
    dat$CO2_flux <- dat$CO2_emis_rate/pars$area
    dat$dCOD_conc_fresh <- pars$conc_fresh$VFA + pars$conc_fresh$Sp + 
        sum(pars$xa_fresh)
    dat$COD_conc_fresh <- pars$conc_fresh$COD
    dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
    dat$ndCOD_conc <- ndCOD_conc <- pars$conc_fresh$COD - dat$dCOD_conc_fresh
    dat$dCOD_conc <- dCOD_conc <- dat$Sp_conc + dat$VFA_conc + 
        rowSums(dat[, paste0(mic_names, "_", "conc")])
    dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
    dat$VS_conc <- pars$COD_conv[["VS"]] * COD_conc
    dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
    dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
    dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
    dat$VS_load_rate <- pars$COD_conv[["VS"]] * dat$COD_load_rate
    dat$COD_load_cum <- cumsum(dat$COD_load_rate * delta_t)
    dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * delta_t)
    dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * delta_t)
    dat$VS_load_cum <- cumsum(dat$VS_load_rate * delta_t)
    dat$CH4_emis_rate_COD <- dat$CH4_emis_rate/dat$COD_load_rate
    dat$CH4_emis_rate_dCOD <- dat$CH4_emis_rate/dat$dCOD_load_rate
    dat$CH4_emis_rate_VS <- dat$CH4_emis_rate/dat$VS_load_rate
    dat$CH4_emis_cum_COD <- dat$CH4_emis_cum/dat$COD_load_cum
    dat$CH4_emis_cum_dCOD <- dat$CH4_emis_cum/dat$dCOD_load_cum
    dat$CH4_emis_cum_VS <- dat$CH4_emis_cum/dat$VS_load_cum
    dat$CO2_emis_rate_COD <- dat$CO2_emis_rate/dat$COD_load_rate
    dat$CO2_emis_rate_dCOD <- dat$CO2_emis_rate/dat$dCOD_load_rate
    dat$CO2_emis_rate_VS <- dat$CO2_emis_rate/dat$VS_load_rate
    dat$CO2_emis_cum_COD <- dat$CO2_emis_cum/dat$COD_load_cum
    dat$CO2_emis_cum_dCOD <- dat$CO2_emis_cum/dat$dCOD_load_cum
    dat$CO2_emis_cum_VS <- dat$CO2_emis_cum/dat$VS_load_cum
    dat$f_COD_CH4_rate <- dat$CH4_emis_rate/pars$COD_conv[["CH4"]]/dat$COD_load_rate
    dat$f_COD_CH4_cum <- dat$COD_conv_cum_meth/dat$COD_load_cum
    dat$f_COD_respir_cum <- dat$COD_conv_cum_respir/dat$COD_load_cum
    dat$f_COD_sr_cum <- dat$COD_conv_cum_sr/dat$COD_load_cum
    names(dat) <- gsub("\\.", "_", names(dat))
    dat_sel <- dat[dat$time > startup, ]
    rCH4_ave <- (dat_sel$CH4_emis_cum[nrow(dat_sel)] - dat_sel$CH4_emis_cum[1])/(dat_sel$time[nrow(dat_sel)] - 
        dat_sel$time[1])
    if (value == "ave") 
        return(rCH4_ave)
    if (value == "ts") 
        return(dat)
    if (value == "tsel") 
        return(dat_sel)
    return(list(pars = pars, ts = dat, tsel = dat_sel, ave = rCH4_ave))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
