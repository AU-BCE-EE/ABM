# =============================================================================
# Test script for the ABM package
# =============================================================================
#
# Runs abm() once for (almost) every argument and argument element, changing
# one thing at a time from the defaults, and checks that
#   * the call succeeds (or fails, where an error is expected),
#   * the output has the expected form (time vector, columns, no NA/Inf),
#   * state variables are non-negative and cumulative emissions never decrease,
#   * conservative components (ash, iNDF) and slurry mass balance,
#   * key results match stored reference values (regression test).
#
# Usage (from the repository root or from the test folder):
#   Rscript test/test_abm.R                    # run all tests on the source code
#   Rscript test/test_abm.R --update-reference # run and overwrite reference values
#   Rscript test/test_abm.R --filter=cover     # only run tests whose name matches 'cover'
#   Rscript test/test_abm.R --installed        # test the installed ABM package instead of the source
#   Rscript test/test_abm.R --verbose          # print warnings for each test
#
# The source code is loaded with pkgload::load_all() (compiles src/ if needed),
# so the package does not have to be re-installed after a change.
#
# When a change to the model is *intended* to change results, check the
# reported differences and then re-run with --update-reference.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
update_ref <- '--update-reference' %in% args
use_installed <- '--installed' %in% args
verbose <- '--verbose' %in% args
filter <- sub('^--filter=', '', grep('^--filter=', args, value = TRUE))
lib_arg <- sub('^--lib=', '', grep('^--lib=', args, value = TRUE))

# Locate test folder and package root
script_dir <- tryCatch({
  f <- sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else stop()
}, error = function(e) {
  if (file.exists('test/test_abm.R')) normalizePath('test') else normalizePath('.')
})
pkg_root <- dirname(script_dir)
ref_file <- file.path(script_dir, 'reference_values.rds')

# Load package
if (use_installed || length(lib_arg) > 0) {
  if (length(lib_arg) > 0) {
    suppressPackageStartupMessages(library(ABM, lib.loc = lib_arg))
  } else {
    suppressPackageStartupMessages(library(ABM))
  }
  cat('Testing installed ABM version', as.character(packageVersion('ABM', lib.loc = if (length(lib_arg)) lib_arg)), '\n')
} else {
  if (!requireNamespace('pkgload', quietly = TRUE)) stop('Package pkgload is needed (or use --installed)')
  suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE, export_all = FALSE))
  cat('Testing ABM source in', pkg_root, '\n')
}
suppressPackageStartupMessages(library(deSolve))

# =============================================================================
# Helpers
# =============================================================================

# Default mng_pars (taken from abm() formals so it is always in sync)
def_mng <- eval(formals(abm)$mng_pars)

# Modify default mng_pars
mng <- function(...) {
  m <- def_mng
  mods <- list(...)
  for (i in names(mods)) m[[i]] <- mods[[i]]
  m
}

# Modify one element of a parameter list
modlist <- function(x, ...) {
  mods <- list(...)
  for (i in names(mods)) x[[i]] <- mods[[i]]
  x
}

# Run abm quietly, capture warnings
run_abm <- function(...) {
  warns <- character()
  res <- withCallingHandlers(
    {
      out <- NULL
      utils::capture.output(out <- suppressMessages(abm(...)))
      out
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart('muffleWarning')
    })
  attr(res, 'warnings') <- unique(warns)
  res
}

check <- function(cond, msg) {
  if (!isTRUE(all(cond))) stop(msg, call. = FALSE)
  invisible(TRUE)
}

state_cols <- c('slurry_mass', 'xa_dead', 'RFd', 'iNDF', 'ash', 'VSd', 'starch', 'CPs', 'CPf',
                'Cfat', 'VFA', 'VFA_C', 'urea', 'TAN', 'sulfate', 'sulfide')
cum_cols <- c('NH3_emis_cum', 'N2O_emis_cum', 'CH4_emis_cum', 'CO2_emis_cum', 'COD_conv_cum',
              'COD_conv_cum_meth', 'COD_conv_cum_respir', 'COD_conv_cum_sr', 'COD_load_cum',
              'C_load_cum', 'N_load_cum', 'slurry_load_cum', 'COD_eff_cum', 'C_eff_cum', 'N_eff_cum')
bal_cols <- as.vector(outer(c('COD', 'C', 'N'), c('in', 'stored', 'removed', 'emitted', 'gap', 'resid', 'rel_resid'),
                            function(e, x) paste0(e, '_bal_', x)))
gap_cols <- c('COD_gap_ferm_cum', 'COD_gap_growth_cum', 'C_gap_ferm_cum', 'C_gap_xa_dead_cum', 'C_gap_uptake_cum',
              'C_gap_respir_cum', 'C_gap_urea_cum', 'N_gap_ferm_cum', 'N_gap_respir_cum')
COD_pools <- c('xa_dead', 'RFd', 'iNDF', 'VSd', 'starch', 'CPs', 'CPf', 'Cfat', 'VFA')
rate_cols <- c('CH4_emis_rate', 'NH3_emis_rate', 'CO2_emis_rate', 'N2O_emis_rate', 'temp_C', 'pH')

# Standard checks on time series output
# grps: microbial group names expected in output
# mass_balance: check slurry mass balance (not valid with washing)
# conserv: check ash/iNDF conservation (only valid with constant conc_fresh)
check_ts <- function(d, days, grps = c('m0', 'm1', 'm2', 'sr1'), times = NULL,
                     mass_balance = TRUE, conserv = TRUE, area = def_mng$area) {

  check(is.data.frame(d), 'Output is not a data frame')
  check(nrow(d) > 1, 'Output has < 2 rows')

  # Time
  check(abs(min(d$time)) < 1e-8, paste('First time is not 0 but', min(d$time)))
  check(abs(max(d$time) - days) < 1e-6, paste('Last time is not days =', days, 'but', max(d$time)))
  check(!is.unsorted(d$time), 'Time is not sorted')
  if (!is.null(times)) {
    check(all(round(unique(d$time), 6) %in% round(sort(unique(c(times, days))), 6)), 'Output contains times not requested')
    check(all(round(sort(unique(c(times, days))), 6) %in% round(d$time, 6)), 'Not all requested times are in output')
  }

  # Columns
  cols <- c(grps, state_cols, cum_cols, bal_cols, gap_cols, rate_cols, paste0(c(grps, state_cols[-1]), '_conc'),
            # Microbial group output is named by group
            as.vector(outer(c('rut_', 'qhat_', 'cum_inhib_', 'H2S_inhib_', 'xa_fresh_'), grps, paste0)))
  miss <- cols[!cols %in% names(d)]
  check(length(miss) == 0, paste('Missing columns:', paste(miss, collapse = ', ')))

  # NA/Inf in key columns
  bad <- cols[sapply(cols, function(x) any(!is.finite(d[[x]])))]
  check(length(bad) == 0, paste('Non-finite values in:', paste(bad, collapse = ', ')))
  # COD/C/N output (except ratios that can be 0/0 legitimately)
  more <- c('COD_conc', 'dCOD_conc', 'VS_conc', 'C_conc', 'N_conc', 'COD_conc_fresh', 'C_conc_fresh', 'N_conc_fresh', 'slurry_depth')
  bad <- more[sapply(more, function(x) any(!is.finite(d[[x]])))]
  check(length(bad) == 0, paste('Non-finite values in:', paste(bad, collapse = ', ')))

  # Non-negative states (allow solver noise)
  for (x in c(grps, state_cols)) {
    tol <- 1e-6 * max(1, max(abs(d[[x]])))
    check(all(d[[x]] >= -tol), paste('Negative values in', x, 'min =', signif(min(d[[x]]), 4)))
  }

  # Cumulative variables never decrease
  for (x in cum_cols) {
    dd <- diff(d[[x]])
    tol <- 1e-6 * max(1, max(abs(d[[x]])))
    check(all(dd >= -tol), paste('Cumulative variable decreases:', x))
  }

  # Emptying: effluent columns present
  check('slurry_mass_eff' %in% names(d), 'slurry_mass_eff column missing')

  # COD, C, and N mass balance: bookkeeping residual (after process stoichiometry gaps) must be ~0
  # Valid for all scenarios (removal is tracked in state variables, so also with washing and output times)
  for (e in c('COD', 'C', 'N')) {
    relerr <- max(abs(d[[paste0(e, '_bal_resid')]])) / max(d[[paste0(e, '_bal_in')]])
    check(relerr < 1e-6, paste0(e, ' mass balance does not close (rel. error ', signif(relerr, 3), ')'))
  }

  # Effluent columns agree with COD removal tracked in COD_eff_cum
  # Removal at the last row may or may not be included in COD_eff_cum (depends on whether row is before or after removal)
  mic_eff <- setdiff(grep('^xa_.+_eff$', names(d), value = TRUE), 'xa_dead_eff')
  COD_eff <- rowSums(d[, c(mic_eff, paste0(COD_pools, '_eff')), drop = FALSE])
  COD_eff[is.na(COD_eff)] <- 0
  a <- sum(COD_eff[-1])
  b <- d$COD_eff_cum[nrow(d)] - d$COD_eff_cum[1]
  tol <- 1e-6 * max(1, a)
  check(abs(a - b) < tol || abs(a - COD_eff[nrow(d)] - b) < tol,
        paste0('Effluent columns do not match COD_eff_cum (', signif(a, 6), ' vs. ', signif(b, 6), ')'))
  # Same for C, using reported effluent C concentration
  C_eff <- d$C_eff_conc * d$slurry_mass_eff
  C_eff[!is.finite(C_eff)] <- 0
  a <- sum(C_eff[-1])
  b <- d$C_eff_cum[nrow(d)] - d$C_eff_cum[1]
  tol <- 1e-6 * max(1, a)
  check(abs(a - b) < tol || abs(a - C_eff[nrow(d)] - b) < tol,
        paste0('C_eff_conc does not match C_eff_cum (', signif(a, 6), ' vs. ', signif(b, 6), ')'))

  # Conservation of ash and iNDF (inert): storage + removed = initial + loaded
  # Effluent is recorded in the row *before* removal, so subtract that row's own effluent
  if (conserv) {
    for (x in c('ash', 'iNDF')) {
      xeff <- d[[paste0(x, '_eff')]]
      xeff[is.na(xeff)] <- 0
      lhs <- d[[x]] + cumsum(xeff) - xeff
      rhs <- d[[x]][1] + d[[paste0('conc_fresh_', x)]][1] * (d$slurry_load_cum - d$slurry_load_cum[1])
      relerr <- max(abs(lhs - rhs)) / max(1, max(abs(rhs)))
      check(relerr < 1e-3, paste0(x, ' not conserved (rel. error ', signif(relerr, 3), ')'))
    }
  }

  # Slurry mass balance: storage + removed = initial + loaded + rain - evap
  if (mass_balance) {
    meff <- d$slurry_mass_eff
    meff[is.na(meff)] <- 0
    lhs <- d$slurry_mass + cumsum(meff) - meff
    pe <- (d$rain - d$evap) * area
    rhs <- d$slurry_mass[1] + (d$slurry_load_cum - d$slurry_load_cum[1]) + c(0, cumsum(diff(d$time) * pe[-1]))
    relerr <- max(abs(lhs - rhs)) / max(1, max(abs(rhs)))
    check(relerr < 1e-3, paste0('Slurry mass balance off (rel. error ', signif(relerr, 3), ')'))
  }

  invisible(TRUE)
}

# Key results stored for regression testing
key_results <- function(d) {
  if (is.list(d) && !is.data.frame(d) && !is.null(d$ts)) d <- d$ts
  if (is.data.frame(d)) {
    n <- nrow(d)
    c(CH4_emis_cum = d$CH4_emis_cum[n], NH3_emis_cum = d$NH3_emis_cum[n], CO2_emis_cum = d$CO2_emis_cum[n],
      N2O_emis_cum = d$N2O_emis_cum[n], COD_load_cum = d$COD_load_cum[n], slurry_mass = d$slurry_mass[n],
      VFA_conc = d$VFA_conc[n], TAN_conc = d$TAN_conc[n], mean_CH4_rate = mean(d$CH4_emis_rate),
      xa_tot = sum(d[n, grepl('_conc$', names(d)) & names(d) %in% paste0(c('m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'p1', 'sr1'), '_conc')]))
  } else {
    d
  }
}

# =============================================================================
# Test registry
# =============================================================================

tests <- list()
add_test <- function(name, fun, regression = TRUE) {
  tests[[name]] <<- list(fun = fun, regression = regression)
}

# Shorthand for standard tests: run abm with arguments, apply check_ts
std_test <- function(name, abm_args, days = abm_args$days, check_args = list(), regression = TRUE,
                     extra = NULL) {
  force(abm_args); force(days); force(check_args); force(extra)
  add_test(name, function() {
    d <- do.call(run_abm, abm_args)
    do.call(check_ts, c(list(d = d, days = days), check_args))
    if (!is.null(extra)) extra(d)
    d
  }, regression = regression)
}

# Base simulation length used in most tests (short, but includes emptying events)
D <- 120

# Variable inputs used in several tests
sm_df <- data.frame(time = c(0, 20, 21, 50, 51, 80, 81, 120),
                    slurry_mass = c(5000, 60000, 5000, 70000, 6000, 80000, 7000, 60000))
temp_df <- data.frame(time = c(0, 30, 60, 90, 120), temp_C = c(10, 15, 22, 18, 12))
pH_df <- data.frame(time = c(0, 40, 80, 120), pH = c(7, 6.5, 7.2, 7))
cf_df <- data.frame(time = c(0, 60, 120), as.data.frame(lapply(man_pars2.0$conc_fresh, function(x) rep(x, 3))))
cf_df$VFA <- c(1.7, 5, 2)
cf_df$TAN <- c(0, 1, 2)
xaf_df <- data.frame(time = c(0, 60, 120), m0 = c(0.06, 0.1, 0.06), m1 = c(0.06, 0.02, 0.06),
                     m2 = c(0.06, 0.06, 0.1), sr1 = c(0.06, 0.06, 0.06))

# ---- Defaults -----------------------------------------------------------------
std_test('defaults_365d', list(days = 365))
std_test('defaults_120d', list(days = D))

# ---- days, delta_t, times -----------------------------------------------------
std_test('days_10', list(days = 10))
std_test('days_730', list(days = 730))
std_test('days_non_integer', list(days = 100.5))
std_test('delta_t_0.5', list(days = D, delta_t = 0.5),
         extra = function(d) check(any(abs(d$time - 10.5) < 1e-8), 'time 10.5 missing with delta_t = 0.5'))
std_test('delta_t_2', list(days = D, delta_t = 2))
std_test('delta_t_0.3_not_multiple', list(days = D, delta_t = 0.3))
std_test('times_regular', list(days = D, times = c(0, 5, 10.5, 42, 60, 100)),
         check_args = list(times = c(0, 5, 10.5, 42, 60, 100), mass_balance = FALSE, conserv = FALSE))
std_test('times_variable_slurry', list(days = D, times = c(0, 10, 25.5, 60, 90),
                                       mng_pars = mng(slurry_mass = sm_df)),
         check_args = list(times = c(0, 10, 25.5, 60, 90), mass_balance = FALSE, conserv = FALSE))

# ---- wthr_pars, evap_pars -----------------------------------------------------
std_test('wthr_pars_no_rain', list(days = D, wthr_pars = modlist(wthr_pars2.0, rain = 0)),
         extra = function(d) check(all(d$rain == 0), 'rain not 0'))
std_test('wthr_pars_high_rain', list(days = D, wthr_pars = modlist(wthr_pars2.0, rain = 5)))
std_test('evap_pars_zero', list(days = D, evap_pars = list(evap = 0)))
std_test('evap_pars_high', list(days = D, evap_pars = list(evap = 4)))

# ---- mng_pars -----------------------------------------------------------------
std_test('mng_slurry_prod_rate_0', list(days = D, mng_pars = mng(slurry_prod_rate = 0)))
std_test('mng_slurry_prod_rate_high', list(days = D, mng_pars = mng(slurry_prod_rate = 9000)))
std_test('mng_slurry_mass_init_0', list(days = D, mng_pars = mng(slurry_mass = 0)))
std_test('mng_slurry_mass_init_high', list(days = D, mng_pars = mng(slurry_mass = 150000)))
std_test('mng_storage_depth', list(days = D, mng_pars = mng(storage_depth = 2)))
std_test('mng_resid_depth', list(days = D, mng_pars = mng(resid_depth = 0.2)))
std_test('mng_floor_area_0', list(days = D, mng_pars = mng(floor_area = 0)),
         extra = function(d) check(all(d$NH3_emis_rate_floor == 0), 'floor NH3 emission with floor_area = 0'))
std_test('mng_area', list(days = D, mng_pars = mng(area = 300)), check_args = list(area = 300))
std_test('mng_empty_int_0_never', list(days = D, mng_pars = mng(empty_int = 0)),
         extra = function(d) check(sum(d$slurry_mass_eff > 0) <= 1, 'emptying with empty_int = 0'))
std_test('mng_empty_int_NA_never', list(days = D, mng_pars = mng(empty_int = NA)))
std_test('mng_empty_int_100', list(days = D, mng_pars = mng(empty_int = 100)))
std_test('mng_empty_int_7', list(days = D, mng_pars = mng(empty_int = 7)))
std_test('mng_temp_C_10', list(days = D, mng_pars = mng(temp_C = 10)),
         extra = function(d) check(all(d$temp_C == 10), 'temp_C not 10'))
std_test('mng_temp_C_35', list(days = D, mng_pars = mng(temp_C = 35)))
std_test('mng_temp_C_below_0', list(days = D, mng_pars = mng(temp_C = -2)))
std_test('mng_temp_C_variable', list(days = D, mng_pars = mng(temp_C = temp_df)),
         extra = function(d) check(abs(d$temp_C[d$time == 30][1] - 15) < 1e-8, 'variable temp_C wrong at t = 30'))
std_test('mng_temp_C_outside_data', list(days = D, mng_pars = mng(temp_C = outside_slurry_temp_NIR)))
std_test('mng_wash_int_14', list(days = D, mng_pars = mng(wash_int = 14, wash_water = 20000, rest_d = 5)),
         check_args = list(mass_balance = FALSE))
std_test('mng_wash_int_rest_0', list(days = D, mng_pars = mng(wash_int = 30, wash_water = 20000, rest_d = 0)),
         check_args = list(mass_balance = FALSE))
std_test('mng_wash_water_0', list(days = D, mng_pars = mng(wash_int = 30, wash_water = 0)))
std_test('mng_wash_rest_d_not_multiple', list(days = D, delta_t = 1, mng_pars = mng(wash_int = 30, wash_water = 20000, rest_d = 2.5)),
         check_args = list(mass_balance = FALSE))
std_test('mng_rest_d_10', list(days = D, mng_pars = mng(wash_int = 40, rest_d = 10)), check_args = list(mass_balance = FALSE))

covers <- c('none', 'tent', 'straw', 'natural_crust', 'clay_pebbles', 'floating_PVC',
            'biocover_porous_sheet', 'corrugated_sheets', 'lid', 'oil', 'peat', 'wood_chips')
for (cv in covers) std_test(paste0('mng_cover_', cv), list(days = D, mng_pars = mng(cover = cv)))
add_test('mng_cover_reduces_NH3', function() {
  d0 <- run_abm(days = D, mng_pars = mng(cover = 'none'))
  d1 <- run_abm(days = D, mng_pars = mng(cover = 'lid'))
  check(tail(d1$NH3_emis_cum, 1) < tail(d0$NH3_emis_cum, 1), 'lid does not reduce NH3 emission')
  check(tail(d1$N2O_emis_cum, 1) > tail(d0$N2O_emis_cum, 1), 'cover does not increase N2O emission')
  d1
})

std_test('mng_resid_enrich_none', list(days = D, mng_pars = mng(resid_enrich = c(OM = 0, mic = 0))))
std_test('mng_resid_enrich_high', list(days = D, mng_pars = mng(resid_enrich = c(OM = 2, mic = 2))))
std_test('mng_resid_enrich_scalar', list(days = D, mng_pars = mng(resid_enrich = 0.9)),
         extra = function(d) {
           d2 <- run_abm(days = D, mng_pars = mng(resid_enrich = c(OM = 0.9, mic = 0.9)))
           check(isTRUE(all.equal(d$CH4_emis_cum, d2$CH4_emis_cum)), 'scalar resid_enrich not used for both OM and mic')
         })
std_test('mng_resid_enrich_mixed', list(days = D, mng_pars = mng(resid_enrich = c(OM = 0, mic = 1.5))))
add_test('mng_resid_enrich_effect', function() {
  d0 <- run_abm(days = D, mng_pars = mng(resid_enrich = c(OM = 0, mic = 0)))
  d1 <- run_abm(days = D, mng_pars = mng(resid_enrich = c(OM = 0, mic = 2)))
  i <- which(d0$time > 43)[1]
  check(d1$m1_conc[i] > d0$m1_conc[i], 'mic enrichment does not increase residual microbes')
  check(abs(d1$RFd_conc[i] - d0$RFd_conc[i]) / d0$RFd_conc[i] < 0.05, 'mic enrichment changed RFd strongly')
  d1
})
std_test('mng_graze', list(days = 365, mng_pars = mng(graze = c(start = 'May', duration = 100, hours_day = 12))),
         extra = function(d) {
           i <- d$time > 125 & d$time < 215
           check(all(abs(d$slurry_prod_rate[i] - 5700 * 0.5) < 1e-6), 'slurry_prod_rate not halved during grazing')
           check(all(abs(d$slurry_prod_rate[d$time < 115] - 5700) < 1e-6), 'slurry_prod_rate reduced outside grazing')
         })
for (s in names(def_mng$scale)) {
  sc <- def_mng$scale
  sc[[s]] <- 1.5
  std_test(paste0('mng_scale_', s), list(days = D, mng_pars = mng(scale = sc)))
}
add_test('mng_scale_order_independent', function() {
  sc <- def_mng$scale
  sc[['qhat_opt']] <- 1.5
  d0 <- run_abm(days = D, mng_pars = mng(scale = sc))
  d1 <- run_abm(days = D, mng_pars = mng(scale = rev(sc)))
  check(isTRUE(all.equal(d0$CH4_emis_cum, d1$CH4_emis_cum)), 'Results depend on order of scale elements')
  d1
}, regression = FALSE)

# ---- Variable slurry mass (abm_variable) ---------------------------------------
std_test('slurry_mass_variable', list(days = D, mng_pars = mng(slurry_mass = sm_df)),
         extra = function(d) {
           check(sum(d$slurry_mass_eff > 0) >= 3, 'Effluent not recorded for variable slurry mass emptying')
           # Emptying at day 20 to 5000 kg, plus 1 day of net rain (rain - evap) on the area
           expect <- 5000 + (d$rain[1] - d$evap[1]) * def_mng$area
           check(abs(d$slurry_mass[d$time == 21][1] - expect) / expect < 0.01, 'slurry_mass after emptying not close to data')
         })
std_test('slurry_mass_variable_beyond_days', list(days = 60, mng_pars = mng(slurry_mass = sm_df)))
std_test('slurry_mass_variable_shorter_than_days', list(days = 150, mng_pars = mng(slurry_mass = sm_df)))
std_test('slurry_mass_variable_start_after_0', list(days = D, mng_pars = mng(slurry_mass = sm_df[-1, ])))
std_test('slurry_mass_variable_with_zero', list(days = D, mng_pars = mng(slurry_mass = transform(sm_df, slurry_mass = c(0, 60000, 0, 70000, 6000, 80000, 7000, 60000)))))
std_test('slurry_mass_variable_wash_water', list(days = D, mng_pars = mng(slurry_mass = cbind(sm_df, wash_water = c(0, 0, 10000, 0, 0, 0, 10000, 0)))),
         check_args = list(mass_balance = FALSE))
for (am in c('early', 'late', 'mid')) {
  std_test(paste0('approx_method_slurry_mass_', am),
           list(days = D, mng_pars = mng(slurry_mass = sm_df), approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = am, conc_fresh = 'constant')),
           check_args = list(mass_balance = am != 'mid'))
}
std_test('approx_method_mid_wash_water', list(days = D, mng_pars = mng(slurry_mass = cbind(sm_df, wash_water = c(0, 0, 10000, 0, 0, 0, 0, 0))),
                                              approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = 'mid', conc_fresh = 'constant')),
         check_args = list(mass_balance = FALSE))

# ---- approx_method for other inputs -------------------------------------------
std_test('approx_method_temp_constant', list(days = D, mng_pars = mng(temp_C = temp_df), approx_method = c(temp = 'constant', pH = 'linear', slurry_mass = 'early', conc_fresh = 'constant')),
         extra = function(d) check(abs(d$temp_C[d$time == 29][1] - 10) < 1e-8, 'constant interpolation of temp_C not used'))
std_test('approx_method_pH_constant', list(days = D, man_pars = modlist(man_pars2.0, pH = pH_df), approx_method = c(temp = 'linear', pH = 'constant', slurry_mass = 'early', conc_fresh = 'constant')),
         extra = function(d) check(abs(d$pH[d$time == 39][1] - 7) < 1e-8, 'constant interpolation of pH not used'))
std_test('approx_method_conc_fresh_linear', list(days = D, man_pars = modlist(man_pars2.0, conc_fresh = cf_df), approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = 'early', conc_fresh = 'linear')),
         check_args = list(conserv = FALSE),
         extra = function(d) check(abs(d$conc_fresh_VFA[d$time == 30][1] - (1.7 + 5) / 2) < 1e-6, 'linear interpolation of conc_fresh not used'))

# ---- man_pars -----------------------------------------------------------------
man_sets <- c('man_pars1.0', 'man_pars2.0', 'man_pars_pig2.0', 'man_pars_pig3.0', 'man_pars_cattle2.0', 'man_pars_cattle3.0')
for (m in man_sets) std_test(paste0('man_pars_', m), list(days = D, man_pars = get(m)))
std_test('man_pars_VSd_set', list(days = D, man_pars = man_pars_VS_pig3.0, grp_pars = grp_pars_VS_pig3.0))
std_test('man_pars_conc_fresh_variable', list(days = D, man_pars = modlist(man_pars2.0, conc_fresh = cf_df)),
         check_args = list(conserv = FALSE),
         extra = function(d) check(abs(d$conc_fresh_VFA[d$time == 70][1] - 5) < 1e-8, 'variable conc_fresh not used'))
std_test('man_pars_pH_6', list(days = D, man_pars = modlist(man_pars2.0, pH = 6)),
         extra = function(d) check(all(d$pH == 6), 'pH not 6'))
std_test('man_pars_pH_variable', list(days = D, man_pars = modlist(man_pars2.0, pH = pH_df)))
std_test('man_pars_pH_calc', list(days = D, man_pars = modlist(man_pars2.0, pH = 'calc')),
         extra = function(d) check(all(d$pH > 4 & d$pH <= 7), 'calculated pH out of range'))
std_test('man_pars_dens', list(days = D, man_pars = modlist(man_pars2.0, dens = 1050)))
std_test('man_pars_conc_fresh_reordered', list(days = D, man_pars = modlist(man_pars2.0, conc_fresh = rev(man_pars2.0$conc_fresh))),
         regression = FALSE,
         extra = function(d) {
           d0 <- run_abm(days = D)
           check(isTRUE(all.equal(d$CH4_emis_cum, d0$CH4_emis_cum)), 'Results depend on order of conc_fresh elements')
         })
add_test('error_VSd_and_other_OM', function() {
  cf <- man_pars2.0$conc_fresh
  cf$VSd <- 50
  r <- tryCatch(run_abm(days = D, man_pars = modlist(man_pars2.0, conc_fresh = cf)), error = function(e) 'error')
  check(identical(r, 'error'), 'No error with both VSd and other OM components > 0')
  NULL
}, regression = FALSE)

# ---- init_pars ----------------------------------------------------------------
std_test('init_pars_conc_init', list(days = D, init_pars = list(conc_init = modlist(man_pars2.0$conc_fresh, VFA = 8, TAN = 2))),
         extra = function(d) check(abs(d$VFA_conc[1] - 8) < 1e-8, 'conc_init not used'))

# ---- grp_pars -----------------------------------------------------------------
grp_sets <- c('grp_pars1.0', 'grp_pars2.0', 'grp_pars_pig2.0', 'grp_pars_pig3.0', 'grp_pars_cattle2.0', 'grp_pars_cattle3.0')
for (g in grp_sets) {
  gp <- get(g)
  std_test(paste0('grp_pars_', g), list(days = D, grp_pars = gp), check_args = list(grps = gp$grps))
}
std_test('grp_pars_subset_groups', list(days = D, grp_pars = modlist(grp_pars2.0, grps = c('m1', 'sr1'))),
         check_args = list(grps = c('m1', 'sr1')))
std_test('grp_pars_single_group', list(days = D, grp_pars = modlist(grp_pars2.0, grps = 'm1')),
         check_args = list(grps = 'm1'))
std_test('grp_pars_no_sr', list(days = D, grp_pars = modlist(grp_pars2.0, grps = c('m0', 'm1', 'm2'))),
         check_args = list(grps = c('m0', 'm1', 'm2')))
std_test('grp_pars_add_m3', list(days = D, grp_pars = modlist(grp_pars2.0, grps = c('m0', 'm1', 'm2', 'm3', 'sr1'),
                                                               xa_fresh = c(default = 0.05))),
         check_args = list(grps = c('m0', 'm1', 'm2', 'm3', 'sr1')))
std_test('grp_pars_xa_fresh_variable', list(days = D, grp_pars = modlist(grp_pars2.0, xa_fresh = xaf_df)),
         check_args = list(conserv = TRUE))
std_test('grp_pars_xa_init', list(days = D, grp_pars = modlist(grp_pars2.0, xa_init = c(all = 0.5))),
         extra = function(d) check(abs(d$m1_conc[1] - 0.5) < 1e-8, 'xa_init not used'))
add_test('grp_pars_decay_rate_per_group', function() {
  # Higher decay rate for m1 only should reduce m1 but (almost) not m2
  d0 <- run_abm(days = D)
  d1 <- run_abm(days = D, grp_pars = modlist(grp_pars2.0, decay_rate = c(default = 0.02, m1 = 0.08)))
  check_ts(d1, D)
  n <- nrow(d1)
  check(d1$m1_conc[n] < 0.9 * d0$m1_conc[n], 'Group-specific decay_rate for m1 has no effect on m1')
  check(abs(d1$m2_conc[n] - d0$m2_conc[n]) / d0$m2_conc[n] < 0.1, 'Group-specific decay_rate for m1 changes m2 strongly')
  d1
})
add_test('wash_enrichment', function() {
  # With washing, residual after the second (wash) emptying is also enriched for particles
  # Enrichment changes composition but not residual slurry mass
  m <- mng(wash_int = 30, wash_water = 20000, rest_d = 0, empty_int = 0)
  d0 <- run_abm(days = 40, mng_pars = modlist(m, resid_enrich = c(OM = 0, mic = 0)))
  d1 <- run_abm(days = 40, mng_pars = modlist(m, resid_enrich = c(OM = 2, mic = 2)))
  i <- which(d0$time == 30)[2]
  check(abs(d1$slurry_mass[i] - d0$slurry_mass[i]) < 1e-6, 'slurry mass after wash differs')
  check(d1$RFd[i] > 1.05 * d0$RFd[i], 'resid_enrich has no effect on residual RFd after wash')
  d1
}, regression = FALSE)
add_test('stoichiometry_balanced', function() {
  # Fermentation and respiration reactions in hard_pars() must be balanced for C, H, O, N, and COD,
  # with substrate COD equal to mol_COD (used in rates_cpp() to convert g COD to mol)
  p <- ABM:::hard_pars(list(pH_LL = c(a = 1), graze_int = 0, graze = c(hours_day = 0), kl = c(NH3 = 1),
                            EF_NH3 = 1, area = 1, EF_N2O = 0))
  atoms <- function(f) {
    m <- regmatches(f, gregexpr('[A-Z][a-z]?[0-9.]*', f))[[1]]
    out <- c(C = 0, H = 0, O = 0, N = 0)
    for (x in m) {
      e <- sub('[0-9.]+$', '', x); n <- sub('^[A-Za-z]+', '', x)
      out[e] <- out[e] + ifelse(n == '', 1, as.numeric(n))
    }
    out
  }
  sub_form <- c(carb = 'C6H10O5', pro = 'C4H6.1O1.2N', lip = 'C51H98O6', OM = 'C15.89H26.3O9.2N')
  # COD (g/mol) of non-substrate species (N as NH3; O2 is negative COD)
  COD_sp <- c(NH3 = 0, H2O = 0, C5H7O2N = 160, C2H4O2 = 64, H2 = 16, CO2 = 0, H2CO3 = 0, O2 = -32)
  # Species in positions 4 to 9 of stoichiometry vectors (positions 1 to 3 are substrate labels)
  sp_ferm <- c('NH3', 'H2O', 'C5H7O2N', 'C2H4O2', 'H2', 'CO2')
  # Respiration: NH4 is released as TAN (NH3), HCO3 counted as CO2 + H2O (H2CO3)
  sp_resp <- c('NH3', 'H2CO3', 'O2', 'H2O', 'C5H7O2N', 'CO2')
  for (s in names(sub_form)) {
    for (type in c('ferm', 'resp')) {
      st <- if (type == 'ferm') p[[s]] else p[[paste0(s, '_resp')]]
      sp <- if (type == 'ferm') sp_ferm else sp_resp
      coef <- st[4:9]
      el <- -atoms(sub_form[[s]])
      for (i in seq_along(sp)) el <- el + coef[[i]] * atoms(sp[i])
      COD <- -p$mol_COD[[s]] + sum(coef * COD_sp[sp])
      # Fermentation is exact (from predFerm()), respiration is rounded
      tol <- if (type == 'ferm') 1e-9 else 1e-3
      check(all(abs(el) < tol * max(1, atoms(sub_form[[s]]))),
            paste0(s, ' ', type, ' not balanced: ', paste(names(el), signif(el, 3), collapse = ', ')))
      check(abs(COD) < tol / 10 * p$mol_COD[[s]], paste0(s, ' ', type, ' COD not balanced (', signif(COD, 4), ' g/mol)'))
      # Fermentation: 2/3 of COD to acetate, 1/3 to H2
      if (type == 'ferm') {
        f_ace <- st[['C2H4O2']] * 64 / (st[['C2H4O2']] * 64 + st[['H2']] * 16)
        check(abs(f_ace - 2/3) < 1e-9, paste0(s, ' fermentation acetate fraction of COD is ', signif(f_ace, 5), ', not 2/3'))
      }
    }
  }
  p[c('carb', 'pro', 'lip', 'OM')]
}, regression = FALSE)
add_test('variable_slurry_mass_rain_exceeds_increase', function() {
  # Given slurry mass increases less than net rain - evaporation in some intervals:
  # slurry addition is set to 0 and rain reduced, so simulated mass follows the given mass (with a warning)
  sm <- data.frame(time = c(0, 50, 51, 120, 121, 250, 251, 365),
                   slurry_mass = c(5000, 70000, 5000, 90000, 6000, 120000, 7000, 80000))
  d <- run_abm(days = 365, mng_pars = mng(slurry_mass = sm), wthr_pars = modlist(wthr_pars2.0, rain = 1.9))
  check_ts(d, 365, mass_balance = TRUE)
  check(any(grepl('increases less than net rain', attr(d, 'warnings'))), 'No warning about rain exceeding slurry mass increase')
  # Simulated mass at end of each filling interval matches given mass
  # (within one day of rain, which falls after emptying with approx_method = "early")
  pe <- (1.9 - d$evap[1]) * def_mng$area
  for (tt in c(50, 120, 250, 365)) {
    sim <- d$slurry_mass[d$time == tt][1]
    given <- sm$slurry_mass[sm$time == tt]
    check(abs(sim - given) <= 1.01 * pe, paste0('Simulated slurry mass at time ', tt, ' (', round(sim), ') differs from given (', given, ')'))
  }
  # With default rain (0) there is no warning
  d0 <- run_abm(days = 365, mng_pars = mng(slurry_mass = sm))
  check(!any(grepl('increases less than net rain', attr(d0, 'warnings'))), 'Unexpected rain warning with default rain = 0')
  d
}, regression = FALSE)
add_test('wash_water_not_at_removal_warns', function() {
  # wash_water is only applied at a removal event: on another row it is ignored, with a warning
  sm_on <- cbind(sm_df, wash_water = c(0, 0, 10000, 0, 0, 0, 0, 0))    # at the emptying (day 21)
  sm_off <- cbind(sm_df, wash_water = c(0, 10000, 0, 0, 0, 0, 0, 0))   # on a filling row (day 20)
  sm_no <- cbind(sm_df, wash_water = 0)
  d_on <- run_abm(days = D, mng_pars = mng(slurry_mass = sm_on, wash_water = 0))
  d_off <- run_abm(days = D, mng_pars = mng(slurry_mass = sm_off, wash_water = 0))
  d_no <- run_abm(days = D, mng_pars = mng(slurry_mass = sm_no, wash_water = 0))
  check(!any(grepl('wash_water value', attr(d_on, 'warnings'))), 'Warning for wash_water at a removal event')
  check(any(grepl('wash_water value.*not applied', attr(d_off, 'warnings'))), 'No warning for wash_water not at a removal event')
  # tolerance: ODE solver noise in slurry mass (~1e-6 kg)
  check(abs(sum(d_on$slurry_mass_eff) - sum(d_no$slurry_mass_eff) - 10000) < 1e-3, 'wash_water at removal not applied')
  check(abs(sum(d_off$slurry_mass_eff) - sum(d_no$slurry_mass_eff)) < 1e-3, 'wash_water not at removal changed effluent')
  w_off <- run_abm(days = D, mng_pars = mng(slurry_mass = sm_off, wash_water = 0), warn = FALSE)
  check(!any(grepl('wash_water value', attr(w_off, 'warnings'))), 'wash_water warning with warn = FALSE')
  d_off
}, regression = FALSE)
add_test('no_emptying_at_end_of_run', function() {
  # Simulation end (100 d) is not an emptying time (42, 84 d), so there must be no effluent at the end
  d <- run_abm(days = 100, mng_pars = mng(empty_int = 42))
  check_ts(d, 100)
  check(d$slurry_mass_eff[nrow(d)] == 0, 'Effluent at end of simulation that is not an emptying time')
  check(sum(d$slurry_mass_eff > 0) == 2, paste('Expected 2 emptying events, found', sum(d$slurry_mass_eff > 0)))
  d
}, regression = FALSE)
add_test('emptying_at_end_when_scheduled', function() {
  # Simulation end (84 d) is a scheduled emptying time, so effluent is reported in the last row
  d <- run_abm(days = 84, mng_pars = mng(empty_int = 42))
  check_ts(d, 84)
  check(d$slurry_mass_eff[nrow(d)] > 0, 'No effluent at end of simulation that is an emptying time')
  d
}, regression = FALSE)
add_test('summ_mass_balance', function() {
  s <- run_abm(days = D, value = 'summ')
  nms <- c('COD_bal_rel_resid', 'C_bal_rel_resid', 'N_bal_rel_resid', 'COD_bal_rel_gap', 'C_bal_rel_gap', 'N_bal_rel_gap')
  check(all(nms %in% names(s)), paste('Missing in summ:', paste(nms[!nms %in% names(s)], collapse = ', ')))
  check(all(is.finite(s[nms])), 'Non-finite mass balance values in summ')
  check(all(abs(s[nms[1:3]]) < 1e-6), 'Mass balance residual in summ not ~0')
  s
}, regression = FALSE)
std_test('grp_pars_decay_rate', list(days = D, grp_pars = modlist(grp_pars2.0, decay_rate = c(all = 0.05))))

# ---- mic_pars, chem_pars, arrh_pars --------------------------------------------
std_test('mic_pars_km_urea', list(days = D, mic_pars = modlist(mic_pars2.0, km_urea = 2)))
std_test('mic_pars_ks_SO4', list(days = D, mic_pars = modlist(mic_pars2.0, ks_SO4 = 0.02)))
std_test('chem_pars_default', list(days = D, chem_pars = chem_pars2.0))
for (a in c('arrh_pars2.0', 'arrh_pars_pig2.0', 'arrh_pars_pig3.0', 'arrh_pars_cattle2.0', 'arrh_pars_cattle3.0')) {
  std_test(paste0('arrh_pars_', a), list(days = D, arrh_pars = get(a)))
}
std_test('arrh_pars_kl_reordered', list(days = D, arrh_pars = modlist(arrh_pars2.0, kl = rev(arrh_pars2.0$kl))),
         regression = FALSE,
         extra = function(d) {
           d0 <- run_abm(days = D)
           check(isTRUE(all.equal(d$NH3_emis_cum, d0$NH3_emis_cum)), 'Results depend on order of kl elements')
         })

# ---- anim_pars ----------------------------------------------------------------
anim_sets <- c('pig_pars2.0', 'pig_pars3.0', 'cattle_pars2.0', 'cattle_pars3.0', 'pig_parsVS2.0', 'pig_parsVS3.0',
               'cattle_parsVS2.0', 'cattle_parsVS3.0', 'digestate_parsVS2.0', 'digestate_parsVS3.0')
for (a in anim_sets) {
  ap <- get(a)
  std_test(paste0('anim_pars_', a), list(days = D, anim_pars = ap), check_args = list(grps = ap$grps))
}

# ---- resp, pH_inhib_overrule ------------------------------------------------------
std_test('resp_FALSE', list(days = D, resp = FALSE),
         extra = function(d) check(all(d$COD_conv_cum_respir < 1e-6), 'respiration with resp = FALSE'))
add_test('resp_TRUE_has_respiration', function() {
  d <- run_abm(days = D, resp = TRUE)
  check(tail(d$COD_conv_cum_respir, 1) > 0, 'no respiration with resp = TRUE')
  d
}, regression = FALSE)
std_test('pH_inhib_overrule_TRUE', list(days = D, pH_inhib_overrule = TRUE))
std_test('pH_inhib_overrule_TRUE_low_pH', list(days = D, pH_inhib_overrule = TRUE, man_pars = modlist(man_pars2.0, pH = 5.5)))

# ---- add_pars, par_key --------------------------------------------------------
std_test('add_pars_simple', list(days = D, add_pars = list(pH = 6.5)),
         extra = function(d) check(all(d$pH == 6.5), 'add_pars pH not used'))
std_test('add_pars_slurry_prod_rate', list(days = D, add_pars = list(slurry_prod_rate = 8000)))
std_test('add_pars_par_element', list(days = D, add_pars = list(qhat_opt.m1 = 2, conc_fresh.VFA = 10)),
         extra = function(d) check(all(abs(d$conc_fresh_VFA - 10) < 1e-8), 'add_pars conc_fresh.VFA not used'))
std_test('add_pars_vector_element', list(days = D, add_pars = list(qhat_opt = c(m1 = 2))))
std_test('add_pars_all_keyword', list(days = D, add_pars = list(yield = c(all = 0.06))))
std_test('add_pars_grps', list(days = D, add_pars = list(grps = c('m1', 'm2'))), check_args = list(grps = c('m1', 'm2')))
std_test('add_pars_temp_C_df', list(days = D, add_pars = list(temp_C = temp_df)))
std_test('add_pars_slurry_mass_df', list(days = D, add_pars = list(slurry_mass = sm_df)))
std_test('par_key_custom', list(days = D, add_pars = list(qhat_opt__m1 = 2), par_key = '__'),
         extra = function(d) {
           d2 <- run_abm(days = D, add_pars = list(qhat_opt.m1 = 2))
           check(isTRUE(all.equal(d$CH4_emis_cum, d2$CH4_emis_cum)), 'par_key = "__" gives different result than "."')
         })
add_test('error_add_pars_bad_name', function() {
  r <- tryCatch(run_abm(days = D, add_pars = list(not_a_par = 1)), error = function(e) 'error')
  check(identical(r, 'error'), 'No error for unknown add_pars name')
  NULL
}, regression = FALSE)

# ---- pars ---------------------------------------------------------------------
add_test('pars_full_list', function() {
  p <- c(wthr_pars2.0, list(evap = 0.5 * ABM:::et(temp_C = 20, pres_kpa = 101, rs = 10)), def_mng, man_pars2.0,
         list(conc_init = man_pars2.0$conc_fresh), grp_pars2.0, mic_pars2.0, chem_pars2.0, arrh_pars2.0,
         list(days = D), resp = TRUE, pH_inhib_overrule = FALSE)
  d <- run_abm(days = D, pars = p)
  check_ts(d, D)
  d0 <- run_abm(days = D)
  check(isTRUE(all.equal(d$CH4_emis_cum, d0$CH4_emis_cum)), 'pars argument gives different result than default')
  d
})

# ---- startup, starting --------------------------------------------------------
std_test('startup_1', list(days = D, startup = 1))
add_test('startup_2_sum', function() {
  s <- run_abm(days = D, startup = 2, value = 'sum')
  check(is.numeric(s) && 'CH4_emis_cum' %in% names(s), 'startup with value = "sum" does not return summary')
  s
})
add_test('startup_changes_start', function() {
  d0 <- run_abm(days = D)
  d1 <- run_abm(days = D, startup = 1)
  check(!isTRUE(all.equal(d0$m1_conc[1], d1$m1_conc[1])), 'startup does not change initial microbial conc')
  d1
}, regression = FALSE)
add_test('startup_rates_calc', function() {
  d <- run_abm(days = D, startup = 1, rates_calc = 'average')
  dt <- c(1, diff(d$time))
  r <- c(0, diff(d$CH4_emis_cum)) / dt
  check(isTRUE(all.equal(d$CH4_emis_rate[dt > 0], r[dt > 0])), 'rates_calc ignored when startup > 0')
  d
}, regression = FALSE)
add_test('startup_anim_pars', function() {
  d0 <- run_abm(days = D, anim_pars = pig_pars3.0)
  d1 <- run_abm(days = D, anim_pars = pig_pars3.0, startup = 1)
  check(!isTRUE(all.equal(d0[nrow(d0) - 5, 'm1_conc'], d1[nrow(d1) - 5, 'm1_conc'])) ||
          !isTRUE(all.equal(d0$m1_conc[1], d1$m1_conc[1])), 'startup has no effect with anim_pars')
  check(!isTRUE(all.equal(d0$m1_conc[1], d1$m1_conc[1])), 'startup initial conditions ignored with anim_pars')
  d1
}, regression = FALSE)
add_test('error_startup_conc_init_add_pars', function() {
  r <- tryCatch(run_abm(days = D, startup = 1, add_pars = list(conc_init.VFA = 3)), error = function(e) 'error')
  check(identical(r, 'error'), 'No error for startup with conc_init in add_pars')
  NULL
}, regression = FALSE)
add_test('starting', function() {
  d0 <- run_abm(days = D)
  d <- run_abm(days = D, starting = d0)
  check_ts(d, D)
  n <- nrow(d0)
  check(abs(d$slurry_mass[1] - d0$slurry_mass[n]) < 1e-6, 'starting slurry_mass not used')
  check(abs(d$VFA[1] - d0$VFA[n]) < 1e-6, 'starting VFA not used')
  check(abs(d$m1_conc[1] - d0$m1_conc[n]) < 1e-6, 'starting microbial conc not used')
  d
})

# ---- value, rates_calc, warn ----------------------------------------------------
add_test('value_sum', function() {
  s <- run_abm(days = D, value = 'sum')
  check(is.numeric(s), 'value = "sum" not numeric')
  need <- c('CH4_emis_cum', 'NH3_emis_cum', 'COD_load', 'f_COD_CH4', 'CH4_emis_rate')
  check(all(need %in% names(s)), 'Missing summary elements')
  check(all(is.finite(s[need])), 'Non-finite summary values')
  s
})
add_test('value_eff', function() {
  e <- run_abm(days = D, value = 'eff')
  check(is.data.frame(e) && nrow(e) >= 2, 'value = "eff" returns no emptying events')
  check(all(e$slurry_mass_eff > 0), 'eff rows without effluent')
  check(!is.unsorted(e$slurry_mass_eff_cum), 'slurry_mass_eff_cum not cumulative')
  e
}, regression = FALSE)
add_test('value_eff_variable_slurry', function() {
  e <- run_abm(days = D, mng_pars = mng(slurry_mass = sm_df), value = 'eff')
  check(is.data.frame(e) && nrow(e) >= 3, 'value = "eff" returns no emptying events for variable slurry mass')
  e
}, regression = FALSE)
add_test('value_all', function() {
  a <- run_abm(days = D, value = 'all')
  check(is.list(a) && all(c('pars', 'ts', 'summ') %in% names(a)), 'value = "all" does not return pars/ts/summ')
  check_ts(a$ts, D)
  a
})
add_test('rates_calc_average', function() {
  d <- run_abm(days = D, rates_calc = 'average')
  check_ts(d, D)
  dt <- c(1, diff(d$time))
  r <- c(0, diff(d$CH4_emis_cum)) / dt
  check(isTRUE(all.equal(d$CH4_emis_rate[dt > 0], r[dt > 0])), 'average CH4 rate not from cumulative emission')
  d
})
std_test('warn_FALSE', list(days = D, warn = FALSE))

# ---- consistency checks -----------------------------------------------------------
add_test('COD_conc_fresh_includes_xa_fresh', function() {
  d <- run_abm(days = D)
  cf <- man_pars2.0$conc_fresh
  dcod <- cf$VFA + cf$xa_dead + cf$RFd + cf$starch + cf$CPs + cf$CPf + cf$Cfat + cf$VSd + sum(grp_pars2.0$xa_fresh)
  check(abs(d$dCOD_conc_fresh[1] - dcod) < 1e-6, paste('dCOD_conc_fresh', d$dCOD_conc_fresh[1], 'expected', dcod))
  NULL
}, regression = FALSE)
add_test('delta_t_does_not_change_results', function() {
  d1 <- run_abm(days = D, delta_t = 1)
  d2 <- run_abm(days = D, delta_t = 0.25)
  r <- abs(tail(d1$CH4_emis_cum, 1) - tail(d2$CH4_emis_cum, 1)) / tail(d1$CH4_emis_cum, 1)
  check(r < 0.01, paste('CH4 emission differs by', signif(100 * r, 3), '% between delta_t = 1 and 0.25'))
  NULL
}, regression = FALSE)
add_test('load_split_in_summary', function() {
  s <- run_abm(days = D, value = 'sum')
  nms <- c('dCOD_load', 'ndCOD_load', 'VS_load', 'Ninorg_load', 'Norg_load', 'CH4_emis_dCOD', 'CH4_emis_VS',
           'CO2_emis_dCOD', 'CO2_emis_VS', 'NH3_emis_Ninorg', 'NH3_emis_Norg', 'N2O_emis_Ninorg', 'N2O_emis_Norg')
  miss <- setdiff(nms, names(s))
  check(length(miss) == 0, paste('Missing from summary:', paste(miss, collapse = ', ')))
  check(abs(s[['dCOD_load']] + s[['ndCOD_load']] - s[['COD_load']]) < 1e-6 * s[['COD_load']], 'dCOD_load + ndCOD_load != COD_load')
  check(abs(s[['Ninorg_load']] + s[['Norg_load']] - s[['N_load']]) < 1e-6 * s[['N_load']], 'Ninorg_load + Norg_load != N_load')
  check(all(s[nms] > 0 | grepl('N2O', nms)), 'Non-positive load or relative emission in summary')
  NULL
}, regression = FALSE)
add_test('emis_slurry_per_kg', function() {
  d <- run_abm(days = D, rates_calc = 'instant')
  s <- run_abm(days = D, value = 'sum')
  n <- nrow(d)
  for (g in c('NH3', 'N2O', 'CH4', 'CO2')) {
    # Rate per kg slurry in storage (g/kg/d)
    check(max(abs(d[[paste0(g, '_emis_rate_slurry')]] * d$slurry_mass - d[[paste0(g, '_emis_rate')]])) < 1e-9 * max(d[[paste0(g, '_emis_rate')]]) + 1e-12,
          paste(g, '_emis_rate_slurry is not per kg slurry'))
    # Cumulative per kg slurry loaded (g/kg), same in time series and summary
    check(abs(d[[paste0(g, '_emis_cum_slurry')]][n] - s[[paste0(g, '_emis_slurry')]]) <= 1e-9 * abs(s[[paste0(g, '_emis_slurry')]]),
          paste(g, '_emis_cum_slurry (ts) differs from', g, '_emis_slurry (summary)'))
  }
  NULL
}, regression = FALSE)
add_test('stoichiometry_gaps_zero', function() {
  # With exact stoichiometry and conversion factors, processes neither create nor destroy COD, C, or N
  for (ap in c('pig_pars2.0', 'pig_parsVS2.0')) {
    s <- run_abm(days = D, anim_pars = get(ap), resp = TRUE, value = 'sum')
    for (e in c('COD', 'C', 'N')) {
      g <- s[[paste0(e, '_bal_rel_gap')]]
      check(abs(g) < 1e-9, paste0(ap, ': ', e, ' balance gap is ', signif(g, 3), ' (stoichiometry or conversion factors not exact)'))
    }
  }
  NULL
}, regression = FALSE)
add_test('compiled_equals_R_function', function() {
  # lsoda() with compiled rates_derivs() must give exactly the same results as with rates_cpp() as R function
  sm <- data.frame(time = c(0, 50, 51, 120), slurry_mass = c(5000, 70000, 5000, 90000))
  for (args in list(list(days = D), list(days = D, mng_pars = mng(slurry_mass = sm), resp = TRUE))) {
    old <- options(ABM.compiled = TRUE); a <- do.call(run_abm, args)
    options(ABM.compiled = FALSE); b <- do.call(run_abm, args); options(old)
    check(identical(names(a), names(b)) && nrow(a) == nrow(b), 'Different output columns or rows')
    num <- names(a)[sapply(a, is.numeric)]
    check(isTRUE(all.equal(a[num], b[num], tolerance = 1e-12, check.attributes = FALSE)), 'Results differ between compiled and R function')
  }
  NULL
}, regression = FALSE)
add_test('N2O_no_emission_from_empty_storage', function() {
  # N2O emission (area * EF_N2O) must fade out when storage is empty, so TAN does not become negative
  sm <- data.frame(time = c(0, 30, 31, 120), slurry_mass = c(1e-10, 1e-10, 2e5, 4e5))
  d <- run_abm(days = D, mng_pars = mng(slurry_mass = sm), add_pars = list(cover = 'tent', evap = 0, rain = 0))
  check(min(d$TAN) > -1e-6, paste('TAN negative:', signif(min(d$TAN), 3)))
  check(max(d$N2O_emis_cum[d$time <= 30]) < 1e-3, 'N2O emitted from empty storage')
  check(d$N2O_emis_cum[nrow(d)] > 0, 'No N2O emission with slurry present')
  check(max(abs(d$N_bal_rel_resid)) < 1e-10, 'N balance does not close')
  NULL
}, regression = FALSE)
add_test('conc_fresh_VFA_C_optional', function() {
  # VFA_C in conc_fresh is optional: default (VFA as acetate) equals explicit VFA / C_VFA
  # A lower VFA_C (VFA with H2 from an upstream ABM run) changes only C, not COD or CH4, and C balance still closes
  cf <- man_pars2.0$conc_fresh
  vc <- cf[['VFA']] / chem_pars2.0$COD_conv[['C_VFA']]
  a <- run_abm(days = D)
  b <- run_abm(days = D, man_pars = modlist(man_pars2.0, conc_fresh = c(cf, VFA_C = vc)))
  d <- run_abm(days = D, add_pars = list(conc_fresh.VFA_C = 0.7 * vc))
  num <- names(a)[sapply(a, is.numeric)]
  check(isTRUE(all.equal(a[num], b[num], tolerance = 1e-12, check.attributes = FALSE)), 'Default VFA_C differs from explicit VFA / C_VFA')
  check(all(abs(a$conc_fresh_VFA_C - vc) < 1e-12), 'Default conc_fresh_VFA_C not VFA / C_VFA')
  check(isTRUE(all.equal(a$CH4_emis_cum, d$CH4_emis_cum, tolerance = 1e-7)), 'VFA_C changed CH4') # only solver noise
  check(isTRUE(all.equal(a$COD_load_cum, d$COD_load_cum, tolerance = 1e-10)), 'VFA_C changed COD load')
  dC <- (a$C_load_cum - d$C_load_cum) - 0.3 * vc * (a$slurry_load_cum - a$slurry_load_cum[1])
  check(max(abs(dC)) < 1e-6 * max(a$C_load_cum), 'C load not reduced by VFA_C difference')
  check(max(abs(d$C_bal_rel_resid), na.rm = TRUE) < 1e-10, 'C balance does not close with given VFA_C')
  check(abs(d$C_conc_fresh[1] - (a$C_conc_fresh[1] - 0.3 * vc)) < 1e-10, 'C_conc_fresh does not use VFA_C')
  NULL
}, regression = FALSE)
add_test('add_pars_unnamed_vector_parameter', function() {
  # Unnamed single value for a parameter with several elements sets all elements (was silently ignored)
  m0 <- mng(resid_enrich = c(OM = 0, mic = 0))
  ref <- run_abm(days = D, mng_pars = m0, value = 'sum')
  a <- run_abm(days = D, add_pars = list(resid_enrich = 0), value = 'sum')
  check(isTRUE(all.equal(a, ref)), 'add_pars resid_enrich = 0 differs from resid_enrich = c(OM = 0, mic = 0)')
  # Unnamed vector of same length replaces by position
  b <- run_abm(days = D, add_pars = list(resid_enrich = c(0, 0)), value = 'sum')
  check(isTRUE(all.equal(b, ref)), 'add_pars resid_enrich = c(0, 0) differs from resid_enrich = c(OM = 0, mic = 0)')
  # Named values still replace only the named elements
  d <- run_abm(days = D, add_pars = list(resid_enrich = c(mic = 0)), value = 'sum')
  e <- run_abm(days = D, mng_pars = mng(resid_enrich = c(OM = 0.9, mic = 0)), value = 'sum')
  check(isTRUE(all.equal(d, e)), 'add_pars with named element does not match')
  # Wrong length or partly named gives an error
  check(inherits(tryCatch(run_abm(days = D, add_pars = list(resid_enrich = c(0, 0, 0))), error = function(e) e), 'error'),
        'No error for unnamed add_pars of wrong length')
  check(inherits(tryCatch(run_abm(days = D, add_pars = list(resid_enrich = c(OM = 0, 0))), error = function(e) e), 'error'),
        'No error for partly named add_pars')
  NULL
}, regression = FALSE)
add_test('ndCOD_conc_equals_iNDF_conc', function() {
  d <- run_abm(days = D)
  check(max(abs(d$ndCOD_conc - d$iNDF_conc)) < 1e-9, 'ndCOD_conc differs from iNDF_conc')
  check(max(abs(d$COD_conc - d$dCOD_conc - d$ndCOD_conc)) < 1e-9, 'COD_conc != dCOD_conc + ndCOD_conc')
  NULL
}, regression = FALSE)

# =============================================================================
# Run tests
# =============================================================================

if (length(filter) > 0) {
  tests <- tests[grepl(filter, names(tests))]
}

ref <- if (file.exists(ref_file)) readRDS(ref_file) else list()
new_ref <- ref

res <- data.frame(test = names(tests), status = NA_character_, seconds = NA_real_, message = '', stringsAsFactors = FALSE)
cat('Running', length(tests), 'tests\n\n')
t_start <- Sys.time()

for (i in seq_along(tests)) {
  nm <- names(tests)[i]
  t0 <- Sys.time()
  out <- tryCatch({
    r <- tests[[i]]$fun()
    list(ok = TRUE, value = r)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  res$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 2)

  if (!out$ok) {
    res$status[i] <- 'FAIL'
    res$message[i] <- gsub('\n', ' ', out$msg)
  } else {
    res$status[i] <- 'PASS'
    # Regression comparison
    if (tests[[i]]$regression && !is.null(out$value)) {
      kr <- key_results(out$value)
      new_ref[[nm]] <- kr
      if (!update_ref) {
        if (is.null(ref[[nm]])) {
          res$status[i] <- 'NEW'
          res$message[i] <- 'No reference values (run with --update-reference)'
        } else {
          common <- intersect(names(kr), names(ref[[nm]]))
          rel <- abs(kr[common] - ref[[nm]][common]) / pmax(abs(ref[[nm]][common]), 1e-8)
          rel[is.na(rel) & is.na(kr[common]) & is.na(ref[[nm]][common])] <- 0
          if (any(is.na(rel) | rel > 1e-4)) {
            res$status[i] <- 'CHANGED'
            bad <- which(is.na(rel) | rel > 1e-4)
            res$message[i] <- paste(sprintf('%s: %.6g -> %.6g', common[bad], ref[[nm]][common][bad], kr[common][bad]), collapse = '; ')
          }
        }
      }
    }
    if (verbose && !is.null(attr(out$value, 'warnings')) && length(attr(out$value, 'warnings'))) {
      res$message[i] <- paste(res$message[i], '[warnings:', paste(attr(out$value, 'warnings'), collapse = ' | '), ']')
    }
  }
  cat(sprintf('%-8s %-45s %6.2fs  %s\n', res$status[i], nm, res$seconds[i], substr(res$message[i], 1, 200)))
}

cat('\n', strrep('=', 70), '\n', sep = '')
cat('Total time:', round(as.numeric(difftime(Sys.time(), t_start, units = 'secs')), 1), 's\n')
print(table(res$status))

if (update_ref) {
  if (length(filter) > 0) {
    cat('\nNote: --update-reference with --filter only updates the filtered tests.\n')
  }
  saveRDS(new_ref, ref_file)
  cat('\nReference values written to', ref_file, '\n')
}

if (any(res$status == 'FAIL')) {
  cat('\nFAILED tests:\n')
  for (i in which(res$status == 'FAIL')) cat(' -', res$test[i], ':', res$message[i], '\n')
}
if (any(res$status == 'CHANGED')) {
  cat('\nCHANGED results (compared to reference values):\n')
  for (i in which(res$status == 'CHANGED')) cat(' -', res$test[i], ':', res$message[i], '\n')
  cat('If these changes are intended, re-run with --update-reference\n')
}

# Non-zero exit status if anything failed (useful for scripts/CI)
if (!interactive() && any(res$status %in% c('FAIL', 'CHANGED'))) quit(status = 1)
