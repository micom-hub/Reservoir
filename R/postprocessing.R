################################################################################
# postprocessing.R
# sourceCpp("postprocessing_fast.cpp")  # call this before sourcing
################################################################################


################################################################################
# add_observation_noise
#applies negbin noise to cases
################################################################################
add_observation_noise <- function(daily_cases) {
  r_val <- rlnorm(1, meanlog = log(50), sdlog = 1.0)

  reported <- numeric(length(daily_cases))
  mask <- daily_cases > 0

  if (any(mask)) {
    mu    <- daily_cases[mask]
    p_nb  <- r_val / (r_val + mu)
    reported[mask]  <- rnbinom(sum(mask), size = r_val, prob = p_nb)
  }

  list(
    reported   = as.integer(reported),
    overdispersion_r = r_val
  )
}


################################################################################
# apply_noise_limited
#
# Lightweight observation model: NegBin overdispersion only, no reporting pipeline; useful for speed
################################################################################
apply_noise_limited <- function(daily_cases_true) {
  noisy <- add_observation_noise(daily_cases_true)
  list(
    daily_cases  = noisy$reported,
    overdispersion_r = noisy$overdispersion_r
  )
}


################################################################################
# add_reporting_effects
#
# Four sequential stages, each active independently with ~80% probability:
#   1. Underreporting ramp - logistic improvement in detection over time
#   2. Weekday effects - lower reporting on weekends
#   3. Lab batch noise - batch accuracy variation, occasional bad batches
#   4. Reporting delays - cases shifted forward with a shrinking delay dist
################################################################################
add_reporting_effects <- function(true_cases, config = NULL) {
  if (is.null(config)) config <- reservoir_config()
  cases    <- as.numeric(true_cases)
  num_days <- length(cases)
  info     <- list()

  apply_reporting <- runif(1) < config$underreport_prob
  apply_weekday   <- runif(1) < config$weekday_prob
  apply_lab       <- runif(1) < config$lab_batch_prob
  apply_delays    <- runif(1) < config$report_delay_prob

  # Shared ramp parameters used by both reporting and delays
  days_to_max    <- sample(30:365, 1)
  outbreak_start <- sample(0:31, 1)
  steepness      <- runif(1, 4, 8)
  midpoint       <- 0.5
  abs_days       <- (seq_len(num_days) - 1) + outbreak_start
  progress       <- pmin(abs_days / days_to_max, 1)

  # Stage 1: underreporting ramp
  if (apply_reporting) {
    init_rate  <- runif(1, 0.05, 0.4)
    final_rate <- runif(1, 0.25, 0.85)
    base_rates <- init_rate + (final_rate - init_rate) *
      (1 / (1 + exp(-steepness * (progress - midpoint))))
    noise_vec  <- rnorm(num_days, 0, 0.05)
    rate_vec   <- pmax(0.1, pmin(0.98, base_rates + noise_vec))
    cases      <- mapply(function(n, p) rbinom(1, as.integer(n), p),
                         as.integer(cases), rate_vec)
    info$reporting <- list(init_rate = init_rate, final_rate = final_rate)
  }

  # Stage 2: weekday effects (vectorized)
  if (apply_weekday) {
    offset       <- sample(0:6, 1)
    weekday_mean <- c(1.2, 1.0, 1.0, 1.0, 0.9, 0.6, 0.4)
    weekday_sd   <- c(0.15, 0.1, 0.1, 0.1, 0.12, 0.2, 0.2)
    wd_idx  <- (seq_len(num_days) - 1 + offset) %% 7 + 1
    factors <- pmax(0.1, pmin(2.0, rnorm(num_days, weekday_mean[wd_idx], weekday_sd[wd_idx])))
    cases   <- cases * factors
    info$weekday <- TRUE
  }

  # Stage 3: lab batch noise
  if (apply_lab) {
    batch_size <- max(10L, rpois(1, 100))
    n_batches  <- num_days %/% batch_size + 1L
    acc        <- rnorm(n_batches, 1.0, 0.08)
    bad        <- runif(n_batches) < 0.005
    acc[bad]   <- runif(sum(bad), 0.7, 0.85)
    p90        <- quantile(cases, 0.9) + 1e-10
    load       <- cases / p90
    penalty    <- exp(-0.2 * pmax(0, load - 0.8))
    for (b in seq_len(n_batches)) {
      st           <- (b - 1) * batch_size + 1
      en           <- min(b * batch_size, num_days)
      cases[st:en] <- cases[st:en] * acc[b] * penalty[st:en]
    }
    info$lab <- list(batch_size = batch_size)
  }

  # stage 4 reporting delay
  if (apply_delays) {
    init_max_delay  <- sample(7:21, 1)
    final_max_delay <- sample(2:7, 1)
    cases      <- apply_reporting_delays_cpp(
      cases, progress,
      init_max_delay, final_max_delay,
      1.0, 4.0, steepness, midpoint
    )
    info$delay <- list(init_max = init_max_delay, final_max = final_max_delay)
  }

  list(
    reported = as.integer(round(pmax(0, cases))),
    info     = info
  )
}


################################################################################
# enforce_no_post_outbreak
#
# zeros out reported cases after the last day with any true cases.
################################################################################
enforce_no_post_outbreak <- function(true_cases, reported) {
  nz <- which(true_cases > 0)
  if (length(nz) == 0) return(reported)
  last_nonzero <- max(nz)
  if (last_nonzero < length(reported))
    reported[(last_nonzero + 1):length(reported)] <- 0L
  reported
}


################################################################################
# add_smaller_reporting_noise
#
# used for hosp and deaths
################################################################################
add_smaller_reporting_noise <- function(values, noise_scale = 1.0) {
  arr       <- as.numeric(values)
  T         <- length(arr)
  apply_rep <- runif(1) < 0.7
  apply_wd  <- runif(1) < 0.4

  if (apply_rep) {
    init_rate  <- pmax(0.05, pmin(0.95, runif(1, 0.5 * noise_scale,
                                               0.8 * noise_scale)))
    final_rate <- runif(1, 0.80, 0.98)
    days_to_max <- sample(30:300, 1)
    x          <- seq_len(T) - 1
    mid        <- days_to_max / 2
    steep      <- runif(1, 3, 6)
    rates      <- init_rate + (final_rate - init_rate) *
      (1 / (1 + exp(-steep * ((x - mid) / days_to_max))))
    rates      <- pmin(1, pmax(0, rates))
    arr        <- mapply(function(n, p) rbinom(1, as.integer(n), p),
                         as.integer(arr), rates)
  }

  if (apply_wd) {
    offset <- sample(0:6, 1)
    for (i in seq_len(T)) {
      wd     <- (i - 1 + offset) %% 7
      factor <- if (wd == 5) runif(1, 0.85, 0.95)
                else if (wd == 6) runif(1, 0.80, 0.90)
                else 1.0
      arr[i] <- arr[i] * factor
    }
  }

  as.integer(round(pmax(0, arr)))
}



################################################################################
# generate_hosp_and_deaths
#
#
# daily_symp_by_age: A x T integer matrix
# hosp_probs, death_probs: length-A numeric vectors
################################################################################
generate_hosp_and_deaths <- function(daily_symp_by_age,
                                      hosp_probs,
                                      death_probs,
                                      hosp_lag_shape,
                                      hosp_lag_scale,
                                      death_lag_shape,
                                      death_lag_scale) {
  A           <- nrow(daily_symp_by_age)
  T           <- ncol(daily_symp_by_age)
  daily_hosp  <- matrix(0L, nrow = A, ncol = T)
  daily_death <- matrix(0L, nrow = A, ncol = T)

  # Population-level convolution
  # really to help for speed
  lag_days <- seq(0, min(T - 1L, 60L))
  h_kern   <- dgamma(lag_days, shape = hosp_lag_shape,  rate = 1 / hosp_lag_scale)
  d_kern   <- dgamma(lag_days, shape = death_lag_shape, rate = 1 / death_lag_scale)
  h_kern   <- h_kern / sum(h_kern)
  d_kern   <- d_kern / sum(d_kern)

  cfr_mult  <- rlnorm(1, 0, 0.4)
  cfr_drift <- runif(1, -0.4, 0.1)
  cfr_time  <- exp(cfr_drift * (seq_len(T) - 1) / max(T - 1, 1))

  for (a in seq_len(A)) {
    symp <- as.numeric(daily_symp_by_age[a, ])

    # Hospitalizations: NegBin daily count convolution
    nb_r_h  <- max(0.5, rlnorm(1, log(2), 1.5))
    mu_h    <- pmax(0, symp * hosp_probs[a])
    n_hosp  <- rnbinom(T, size = nb_r_h, mu = mu_h)
    hosp_sig <- as.numeric(convolve(n_hosp, rev(h_kern), type = "open"))[seq_len(T)]
    daily_hosp[a, ] <- pmax(0L, as.integer(round(hosp_sig)))

    # Deaths: fraction of hospitalizations convolution
    eff_cfr   <- pmin(0.95, death_probs[a] * cfr_mult * cfr_time)
    nb_r_d    <- max(0.5, rlnorm(1, log(1), 1.5))
    mu_d      <- pmax(0, n_hosp * eff_cfr)
    n_die     <- rnbinom(T, size = nb_r_d, mu = mu_d)
    death_sig <- as.numeric(convolve(n_die, rev(d_kern), type = "open"))[seq_len(T)]
    daily_death[a, ] <- pmax(0L, as.integer(round(death_sig)))
  }

  list(hosp = daily_hosp, death = daily_death)
}



################################################################################
# compute_hospital_occupancy
#
# For each admission draws a LOS from Gamma(los_shape, los_scale) and
# increments occupancy for each day of the stay.
################################################################################
compute_hospital_occupancy <- function(daily_hosp, avg_los = 7L, ...) {
  hospital_occupancy_cpp(as.integer(daily_hosp), as.integer(avg_los))
}


################################################################################
# generate_wastewater
#
# Convolves incidence with a gamma shedding kernel
################################################################################
generate_wastewater <- function(daily_cases,
                                 noise_sigma   = 0.2,
                                 smooth_window = 3L,
                                 pmmov_mean    = 1e6,
                                 pmmov_sigma   = 0.15) {
  T               <- length(daily_cases)
  mean_shed_days  <- runif(1, 3, 7)
  shed_dispersion <- 2.0

  x_kern <- seq(0, 29)
  kernel <- dgamma(x_kern, shape = shed_dispersion,
                   rate  = shed_dispersion / mean_shed_days)
  kernel <- kernel / sum(kernel)

  pathogen <- convolve(as.numeric(daily_cases), rev(kernel),
                       type = "open")[seq_len(T)]

  if (smooth_window > 1L) {
    k        <- rep(1, smooth_window) / smooth_window
    pathogen <- as.numeric(stats::filter(pathogen, k, sides = 2))
    pathogen[is.na(pathogen)] <- 0
  }

  pathogen <- pmax(0, pathogen * rlnorm(T, 0, noise_sigma))
  pmmov    <- pmmov_mean * rlnorm(T, 0, pmmov_sigma)

  list(
    raw        = pathogen,
    normalized = pathogen / pmax(pmmov, 1e-9),
    info       = list(
      mean_shed_days = mean_shed_days,
      noise_sigma    = noise_sigma,
      pmmov_mean     = pmmov_mean
    )
  )
}


################################################################################
# generate_syndromic
# Synthetic syndromic surveillance (e.g. ILI consultations)
################################################################################
generate_syndromic <- function(daily_symp,
                                consult_prob  = 0.1,
                                noise_sigma   = 0.2,
                                smooth_window = 3L) {
  T        <- length(daily_symp)
  expected <- daily_symp * consult_prob
  noisy    <- pmax(0, expected * rlnorm(T, 0, noise_sigma))

  if (smooth_window > 1L) {
    k     <- rep(1, smooth_window) / smooth_window
    noisy <- as.numeric(stats::filter(noisy, k, sides = 2))
    noisy[is.na(noisy)] <- 0
  }

  list(
    signal = as.integer(round(noisy)),
    info   = list(consult_prob = consult_prob, noise_sigma = noise_sigma)
  )
}


################################################################################
# run_postprocessing
#
# Wraps all the signals together. Takes C++ simulator output and runs the full pipeline:
#  1. Observation noise (NegBin)
#   2. Reporting effects (ramp, weekday, lab, delay)
#   3. Post-outbreak zeroing
#   4. Symptomatic extraction
#   5. Hospitalizations and deaths by age
#   6. Hospital occupancy
#   7. Reported hosp and deaths (smaller noise)
#   8. Wastewater signal
#   9. Syndromic surveillance
#
################################################################################
run_postprocessing <- function(result,
                                params,
                                P,
                                A,
                                num_days,
                                age_hosp_probs  = NULL,
                                age_death_probs = NULL,
                                skip_occupancy  = FALSE,
                                config          = NULL) {
  if (is.null(config)) config <- reservoir_config()
  block <- P * A

  # aggregate daily
  daily_true <- rowSums(result$daily_new_cases_block)

  # 1. Observation noise
  noisy<- add_observation_noise(daily_true)
  daily_noisy <- noisy$reported

  # 2. Reporting effects
  rep_out   <- add_reporting_effects(daily_noisy, config)
  daily_rep <- enforce_no_post_outbreak(daily_true, rep_out$reported)

  # 3. Build a x t matrix
  # daily_new_cases_block is num_days x (P*A), columns ordered p*A + a (0-based)
  if (A > 1) {
    # sum acros patches
    cases_mat    <- result$daily_new_cases_block  # num_days x PA
    daily_by_age <- matrix(0, nrow = A, ncol = num_days)
    for (a in seq_len(A)) {
      # Columns belonging to age a 
      age_cols <- seq(a, block, by = A)
      daily_by_age[a, ] <- rowSums(cases_mat[, age_cols, drop = FALSE])
    }
  } else {
    daily_by_age <- matrix(daily_true, nrow = 1, ncol = num_days)
  }

  daily_symp_by_age <- daily_by_age * (1 - params$p_asym)

  # 4. Hospitalizations and deaths
  if (is.null(age_hosp_probs))
    age_hosp_probs  <- seq(0.003, 0.25, length.out = A)
  if (is.null(age_death_probs))
    age_death_probs <- seq(0.001, 0.15, length.out = A)

  # Ensure vectors for A=1 scalar inputs
  age_hosp_probs  <- rep_len(age_hosp_probs,  A)
  age_death_probs <- rep_len(age_death_probs, A)

  if (!isTRUE(config$compute_hosp_deaths)) {
    hd <- list(hosp  = matrix(0L, nrow = A, ncol = num_days),
               death = matrix(0L, nrow = A, ncol = num_days))
  } else hd <- generate_hosp_and_deaths(
    daily_symp_by_age,
    hosp_probs  = age_hosp_probs,
    death_probs = age_death_probs,
    hosp_lag_shape  = runif(1, 3, 6),
    hosp_lag_scale  = runif(1, 1, 3),
    death_lag_shape = runif(1, 1, 4),
    death_lag_scale = runif(1, 2, 8)
  )

  daily_hosp_true  <- colSums(hd$hosp)
  daily_death_true <- colSums(hd$death)

  # 5. Hospital occupancy
  occupancy <- if (skip_occupancy || !isTRUE(config$compute_hosp_deaths))
    integer(num_days) else compute_hospital_occupancy(daily_hosp_true)

  # 6. Reported hosp and deaths
  daily_hosp_rep  <- add_smaller_reporting_noise(daily_hosp_true,  1.0)
  daily_death_rep<- add_smaller_reporting_noise(daily_death_true, 0.75)

  # 7.Wastewater
  ww <- if (isTRUE(config$compute_wastewater)) generate_wastewater(daily_true)
        else list(raw = numeric(num_days), normalized = numeric(num_days), info = list())

  # 8. Syndromic
  synd <- if (isTRUE(config$compute_syndromic))
    generate_syndromic(colSums(daily_symp_by_age))
  else list(signal = integer(num_days), info = list())

  list(
    daily_cases_true      = daily_true,
    daily_cases_reported  = daily_rep,
    daily_hosp_true       = daily_hosp_true,
    daily_hosp_reported   = daily_hosp_rep,
    daily_death_true      = daily_death_true,
    daily_death_reported  = daily_death_rep,
    hosp_occupancy        = occupancy,
    wastewater_raw        = ww$raw,
    wastewater_normalized = ww$normalized,
    syndromic_signal      = synd$signal,
    noise_info = list(
      overdispersion_r = noisy$overdispersion_r,
      reporting        = rep_out$info,
      wastewater       = ww$info
    )
  )
}
