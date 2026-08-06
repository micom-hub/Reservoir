###############################################################################
# setup_sim.R
#
# Builds the cfg list consumed by the C++ engines.
#
###############################################################################


# (.HUMAN_COMPS / .N_HUMAN_BLOCKS are defined in aaa_constants.R)


#' Build the simulation configuration object
#'
#' @param N Total population size.
#' @param P Number of patches. @param A Number of age groups.
#' @param tf_days Simulation horizon in days.
#' @param beta Transmission rate per contact per day.
#' @param infectious_days,latent_days,immunity_days Mean  times. Use
#'   `Inf` for `immunity_days` to make R absorbing (permanent immunity), and
#'   `Inf` for `infectious_days` to disable recovery (SI models).
#' @param p_asym Fraction of infections that are asymptomatic (0 disables Ia).
#' @param asymp_trans Relative infectiousness of asymptomatic cases.
#' @param kappa Rate at which asymptomatic cases become symptomatic.
#' @param n_waning_stages 1, 2, or 3 sequential R compartments.
#' @param no_latent STRUCTURAL. TRUE routes infections directly into Is/Ia,
#'   bypassing L entirely (SIR / SIAR family).
#' @param sis_mode STRUCTURAL. TRUE retargets recovery to S instead of R
#'   (SIS / SLIS family).
#' @param prep_start_rate Daily rate at which susceptibles enter Sprep. This
#'   is the vaccination/PrEP campaign rate — nonzero values are required for
#'   Sprep to populate at all.
#' @param reactive TRUE enables the in-engine threshold-triggered controller.
#'   In reactive mode the intervention driver vectors are supplied at their
#'   ACTIVE values for all days and the engine gates them on/off.
#'
#' @export
setup_sim <- function(

  # Population
  N = 10000L, P = 1L, A = 1L, tf_days = 365L,
  pop_fracs = NULL, patch_fracs = NULL,

  # Epid params
  beta = 0.3,
  infectious_days = 7,
  latent_days     = 5,
  immunity_days   = 365,
  p_asym          = 0.3,
  asymp_trans     = 0.5,
  kappa           = 0.1,
  n_waning_stages = 1L,

  # Structural stuff
  no_latent = FALSE,
  sis_mode  = FALSE,

  # Seeds
  n_seed = 10L, seed_comp = "Is",

  # Demographics
  birth_rate = 0.0, death_rate = 0.0, importation_rate = 0.0,
  spost_waning_rate = 1/14,

  # Waves
  wave_change_days = NULL, wave_change_betas = NULL,

  # Seasonality
  use_seasonality = FALSE, n_harmonics = 1L,
  harmonic_amps = 0.3, harmonic_offsets = 0, harmonic_periods = 365,
  annual_jitter_sd = 0, daily_noise_sd = 0, daily_noise_clamp = c(0.5, 2.0),

  # Contact reduc
  contact_mult = 1.0,

  # Interventions
  iso_rate  = 0.0, iso_on_day  = NULL, iso_off_day  = NULL, iso_active  = 0.3,
  quar_rate = 0.0, quar_on_day = NULL, quar_off_day = NULL, quar_active = 0.2,
  pep_rate  = 0.0, pep_on_day  = NULL, pep_off_day  = NULL, pep_active  = 0.1,
  prep_eff  = 0.0, prep_on_day = NULL, prep_off_day = NULL, prep_active = 0.7,
  leave_quar_rate = 0.1, leave_quar_prep_rate = 0.1,
  prep_start_rate = 0.0, prep_stop_rate = 0.0,

  # Intervention reactiion
  reactive = FALSE,
  on_threshold = 0.0, off_threshold = 0.0, trigger_delay = 0L,
  reactive_contact_mult = 1.0,

  # Matrics
  C = NULL, Mp = NULL,

  # Tau-leaping
  epsilon = 0.03, Ncritical = 10L, exactThreshold = 10.0, maxtau = Inf
) {

  # Validation
  if (N <= 0)       stop("N must be > 0")
  if (P <= 0)       stop("P must be > 0")
  if (A <= 0)       stop("A must be > 0")
  if (tf_days <= 0) stop("tf_days must be > 0")
  if (!n_waning_stages %in% 1:3) stop("n_waning_stages must be 1, 2, or 3")

  block <- P * A

  # Initial states
  comp_offsets <- setNames(seq_along(.HUMAN_COMPS) - 1L, .HUMAN_COMPS)
  if (!seed_comp %in% names(comp_offsets))
    stop("seed_comp must be one of: ", paste(.HUMAN_COMPS, collapse = ", "))
  if (no_latent && seed_comp %in% c("L", "QL"))
    stop("seed_comp = '", seed_comp, "' is invalid when no_latent = TRUE ",
         "(that compartment is unreachable). Use 'Is' or 'Ia'.")

  # Distribute the population across (patch, age) cells.
  if (is.null(pop_fracs))   pop_fracs   <- rep(1/A, A)
  if (is.null(patch_fracs)) patch_fracs <- rep(1/P, P)
  if (length(pop_fracs) != A)
    stop("`pop_fracs` must have length A (", A, ")")
  if (length(patch_fracs) != P)
    stop("`patch_fracs` must have length P (", P, ")")
  pop_fracs   <- pop_fracs   / sum(pop_fracs)
  patch_fracs <- patch_fracs / sum(patch_fracs)

  cell_frac <- rep(patch_fracs, each = A) * rep(pop_fracs, times = P)


  raw    <- N * cell_frac
  n_cell <- floor(raw)
  short  <- N - sum(n_cell)
  if (short > 0) {
    ord <- order(raw - n_cell, decreasing = TRUE)
    n_cell[ord[seq_len(short)]] <- n_cell[ord[seq_len(short)]] + 1
  }

  x0 <- numeric(.N_HUMAN_BLOCKS * block)
  x0[1:block] <- n_cell

  # Seed into cell (0, 0)- a point introduction that then spreads through
  # the mixing and contact matrices.
  n_seed <- min(n_seed, x0[1])
  x0[1]  <- x0[1] - n_seed
  x0[comp_offsets[seed_comp] * block + 1L] <- n_seed

  # wave
  if (!is.null(wave_change_days)) {
    if (is.null(wave_change_betas))
      stop("wave_change_betas must be supplied when wave_change_days is set")
    if (length(wave_change_betas) != length(wave_change_days))
      stop("wave_change_days and wave_change_betas must be the same length")

    all_days  <- c(0L, as.integer(wave_change_days))
    all_betas <- c(beta, as.numeric(wave_change_betas))
    beta_vec  <- numeric(tf_days)
    for (day in seq_len(tf_days)) {
      seg <- sum(all_days < day)
      beta_vec[day] <- all_betas[seg]
    }
  } else {
    beta_vec <- rep(beta, tf_days)
  }

  # Seasonality
  if (use_seasonality && n_harmonics > 0) {
    h_amps    <- rep_len(as.numeric(harmonic_amps), n_harmonics)
    h_offsets <- rep_len(as.numeric(harmonic_offsets), n_harmonics)
    h_periods <- rep_len(as.numeric(harmonic_periods), n_harmonics)

    n_years     <- tf_days %/% 365L + 2L
    year_jitter <- if (annual_jitter_sd > 0)
      rnorm(n_years, 0, annual_jitter_sd) else rep(0, n_years)
    daily_noise <- if (daily_noise_sd > 0)
      pmax(daily_noise_clamp[1],
           pmin(daily_noise_clamp[2], rnorm(tf_days, 1, daily_noise_sd)))
      else rep(1.0, tf_days)

    seasonal_vec <- numeric(tf_days)
    for (d in seq_len(tf_days)) {
      day      <- d - 1L
      year_idx <- day %/% 365L + 1L
      jitter   <- year_jitter[min(year_idx, n_years)]
      val      <- 0.0
      for (h in seq_len(n_harmonics))
        val <- val + h_amps[h] *
          cos(2 * pi * (day - h_offsets[h] - jitter) / h_periods[h])
      seasonal_vec[d] <- pmax(0, 1 + val) * daily_noise[d]
    }
  } else {
    seasonal_vec <- if (daily_noise_sd > 0)
      pmax(daily_noise_clamp[1],
           pmin(daily_noise_clamp[2], rnorm(tf_days, 1, daily_noise_sd)))
      else rep(1.0, tf_days)
  }

  # Contact reduction (school/work closure) folds straight into beta.
  beta_vec <- beta_vec * seasonal_vec * contact_mult

  # Driver helpers
  make_driver <- function(baseline, on_day, off_day, active) {
    if (reactive) return(rep(as.numeric(active), tf_days))
    vec <- rep(as.numeric(baseline), tf_days)
    if (!is.null(on_day)) {
      on_day  <- as.integer(on_day)
      off_day <- if (!is.null(off_day)) as.integer(off_day) else tf_days
      if (on_day < 1 || on_day > tf_days)
        stop("on_day must be between 1 and tf_days (", tf_days, ")")
      if (off_day < on_day) stop("off_day must be >= on_day")
      off_day <- min(off_day, tf_days)
      vec[on_day:off_day] <- as.numeric(active)
    }
    vec
  }

  expand <- function(x, nm) {
    if (length(x) == 1L) return(rep(as.numeric(x), tf_days))
    if (length(x) != tf_days)
      stop(nm, " must be length 1 or tf_days (", tf_days, ")")
    as.numeric(x)
  }

  # Matrices
  if (is.null(C))  C  <- matrix(1.0, A, A)
  if (is.null(Mp)) Mp <- matrix(1.0 / P, P, P)
  if (!is.matrix(C)  || nrow(C)  != A || ncol(C)  != A)
    stop("C must be an A x A matrix (A = ", A, ")")
  if (!is.matrix(Mp) || nrow(Mp) != P || ncol(Mp) != P)
    stop("Mp must be a P x P matrix (P = ", P, ")")

  # Rates
  # immunity_days = Inf  -> omega = 0 -> R absorbing (permanent immunity)
  # infectious_days= Inf -> gamma = 0 -> no recovery (SI)
  omega <- if (is.finite(immunity_days))   1 / max(immunity_days, 1e-8)   else 0.0
  gamma <- if (is.finite(infectious_days)) 1 / max(infectious_days, 1e-8) else 0.0
  sigma <- 1 / max(latent_days, 1e-8)

  list(
    P = P, A = A, tf_days = tf_days,
    n_waning_stages = as.integer(n_waning_stages),

    # structural flags consumed by the C++ engines
    no_latent = as.integer(isTRUE(no_latent)),
    sis_mode  = as.integer(isTRUE(sis_mode)),

    C = C, Mp = Mp, x0 = x0,

    beta_eff_daily        = beta_vec,
    importation_daily     = rep(as.numeric(importation_rate), tf_days),
    iso_rate_daily        = make_driver(iso_rate,  iso_on_day,  iso_off_day,  iso_active),
    quar_rate_daily       = make_driver(quar_rate, quar_on_day, quar_off_day, quar_active),
    leave_quar_daily      = expand(leave_quar_rate,      "leave_quar_rate"),
    leave_quar_prep_daily = expand(leave_quar_prep_rate, "leave_quar_prep_rate"),
    prep_start_daily      = expand(prep_start_rate,      "prep_start_rate"),
    prep_stop_daily       = expand(prep_stop_rate,       "prep_stop_rate"),
    pep_rate_daily        = make_driver(pep_rate,  pep_on_day,  pep_off_day,  pep_active),
    prep_eff_daily        = make_driver(prep_eff,  prep_on_day, prep_off_day, prep_active),

    # reactive controller
    reactive_mode         = as.integer(isTRUE(reactive)),
    on_threshold          = as.numeric(on_threshold),
    off_threshold         = as.numeric(off_threshold),
    trigger_delay         = as.integer(trigger_delay),
    reactive_contact_mult = as.numeric(reactive_contact_mult),

    params = list(
      sigma = sigma, gamma = gamma, omega = omega,
      p_asym = p_asym, asymp_trans = asymp_trans, kappa = kappa,
      birth_rate = birth_rate, death_rate = death_rate,
      spost_waning_rate = spost_waning_rate
    ),

    epsilon = epsilon, Ncritical = Ncritical,
    exactThreshold = exactThreshold, maxtau = maxtau,

    N = N,

    # layout metadata for downstream tooling
    comp_names     = .HUMAN_COMPS,
    n_human_blocks = .N_HUMAN_BLOCKS
  )
}


################################################################################
# set_initial_conditions
#
# Populates cfg$x0 for novel vs non-novel pathogens.
#
# Novel:     fraction-based seeding across L, Is, Ia (and Iso/QL if an
#            intervention is active) proportional to mean sojourn time.
#            Never seeds R, Sprep, Spost, or cumul.
# Non-novel: immune_frac supplied directly, or derived from the endemic
#            equilibrium 1 - 1/R0. Remaining susceptibles seeded as above.
#
# Respects cfg$no_latent: when the latent stage is bypassed, seeding weight
# that would have gone to L / QL is redistributed to Is / Ia.
################################################################################
#' Populate the initial state vector
#'
#' @export
set_initial_conditions <- function(cfg,
                                   novel            = TRUE,
                                   immune_frac      = NULL,
                                   R0               = NULL,
                                   seed_frac        = NULL,
                                   seed_frac_lo     = 1e-6,
                                   seed_frac_hi     = 5e-4,
                                   infectious_days  = 5.0,
                                   latent_days      = 5.0,
                                   p_asym           = 0.30,
                                   has_intervention = FALSE,
                                   waning_stages    = 1L) {

  P <- cfg$P; A <- cfg$A; N <- cfg$N
  block <- P * A
  x0    <- cfg$x0
  idx   <- function(blk, p = 0, a = 0) blk * block + p * A + a + 1L
  no_latent <- isTRUE(cfg$no_latent == 1L)

  #  Pre-existing immunity
  if (novel) {
    immune_frac <- 0.0
  } else if (is.null(immune_frac)) {
    immune_frac <- if (is.null(R0) || is.na(R0) || R0 <= 1) 0.0
                   else rnorm(1, 1 - 1/R0, 0.05)
  }
  immune_frac <- pmax(0, pmin(as.numeric(immune_frac), 0.97))

  # Under SIS there is no R compartment to hold immunity.
  if (isTRUE(cfg$sis_mode == 1L)) immune_frac <- 0.0

  n_R_blocks <- min(waning_stages, 3L)
  if (immune_frac > 0 && n_R_blocks >= 1) {
    n_immune <- round(N * immune_frac)
    per_R    <- round(n_immune / n_R_blocks)
    for (r_blk in 9:(9 + n_R_blocks - 1)) {
      per_cell <- round(per_R / block)
      for (p in 0:(P - 1)) for (a in 0:(A - 1)) {
        s_col <- idx(0, p, a); r_col <- idx(r_blk, p, a)
        transfer <- min(per_cell, x0[s_col])
        x0[s_col] <- x0[s_col] - transfer
        x0[r_col] <- x0[r_col] + transfer
      }
    }
  }

  # ---- Seed active infections (fraction-based, scales with N) ----
  if (is.null(seed_frac))
    seed_frac <- exp(runif(1, log(seed_frac_lo), log(seed_frac_hi)))
  n_total_seed <- max(1L, round(N * seed_frac))

  if (no_latent) {
    # No L / QL to seed into: split entirely between Is and Ia.
    w_Is  <- (1 - p_asym) * infectious_days
    w_Ia  <- p_asym       * infectious_days
    w_Iso <- if (has_intervention) w_Is * 0.10 else 0
    weights <- c(Is = w_Is, Ia = w_Ia, Iso = w_Iso)
    blk_map <- c(Is = 6L, Ia = 7L, Iso = 8L)
  } else {
    w_L   <- latent_days
    w_Is  <- (1 - p_asym) * infectious_days
    w_Ia  <- p_asym       * infectious_days
    w_Iso <- if (has_intervention) w_Is * 0.10 else 0
    w_QL  <- if (has_intervention) w_L  * 0.08 else 0
    weights <- c(L = w_L, Is = w_Is, Ia = w_Ia, Iso = w_Iso, QL = w_QL)
    blk_map <- c(L = 4L, Is = 6L, Ia = 7L, Iso = 8L, QL = 5L)
  }

  if (sum(weights) <= 0) weights[1] <- 1
  weights <- weights / sum(weights)

  for (nm in names(weights)) {
    n_this <- round(n_total_seed * weights[[nm]])
    if (n_this <= 0) next
    s_col <- idx(0L, 0L, 0L)
    c_col <- idx(blk_map[[nm]], 0L, 0L)
    transfer <- min(n_this, x0[s_col])
    x0[s_col] <- x0[s_col] - transfer
    x0[c_col] <- x0[c_col] + transfer
  }

  cfg$x0          <- x0
  cfg$immune_frac <- immune_frac
  cfg$seed_frac   <- seed_frac
  cfg$n_seeded    <- n_total_seed
  cfg
}
