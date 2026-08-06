################################################################################
# config.R
#
################################################################################


################################################################################
# .as_prob
################################################################################-
.as_prob <- function(x, nm) {
  if (is.logical(x)) return(if (isTRUE(x)) 1.0 else 0.0)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1)
    stop("`", nm, "` must be TRUE/FALSE or a probability in [0, 1]")
  as.numeric(x)
}

################################################################################
# .draw_flag
################################################################################
.draw_flag <- function(p) {
  if (p >= 1) return(TRUE)
  if (p <= 0) return(FALSE)
  runif(1) < p
}


################################################################################
# model_structure
#
# Declares which compartments and transitions exist. This is the object that decides
# the compartmental skeleton (SIR / SEIR / SIS / SEIAR / ...).
#
# Each flag accepts TRUE / FALSE / probability.
################################################################################
#' Declare the compartmental model structure
#'
#' @param has_latent Latent (exposed) stage present. STRUCTURAL: when FALSE,
#'   infections deposit directly into Is/Ia and L/QL stay empty.
#' @param has_asymptomatic Asymptomatic infectious class. Parametric (p_asym = 0).
#' @param has_recovery Recovery from infection. FALSE gives SI / SLI (chronic
#'   infection); requires `beta` supplied directly since R0/IP inversion diverges.
#' @param immunity One of "waning" (R -> S at rate omega), "permanent"
#'   (omega = 0, R absorbing), or "none" (STRUCTURAL SIS: recovery returns
#'   directly to S and R stays empty).
#' @param n_waning_stages 1, 2, or 3 sequential R compartments.
#' @param has_isolation Isolation of infectious individuals into Iso.
#' @param has_quarantine Quarantine of susceptible (Qs) and latent (QL).
#' @param has_pep Post-exposure prophylaxis. Requires `has_latent` (acts on L/QL).
#' @param has_vaccination Leaky vaccination via the Sprep compartment.
#' @param has_demography Background births and deaths (population turnover).
#' @param has_importation External case importation.
#'
#' @export
model_structure <- function(
    has_latent = TRUE,
    has_asymptomatic = TRUE,
    has_recovery = TRUE,
    immunity = c("waning", "permanent", "none"),
    n_waning_stages  = 1,
    has_isolation = TRUE,
    has_quarantine = FALSE,
    has_pep = FALSE,
    has_vaccination = FALSE,
    has_demography = TRUE,
    has_importation = TRUE
) {
  immunity <- match.arg(immunity)

  if (!n_waning_stages %in% 1:3)
    stop("`n_waning_stages` must be 1, 2, or 3")

  s <- list(
    has_latent = .as_prob(has_latent,       "has_latent"),
    has_asymptomatic = .as_prob(has_asymptomatic, "has_asymptomatic"),
    has_recovery = .as_prob(has_recovery,     "has_recovery"),
    immunity = immunity,
    n_waning_stages  = as.integer(n_waning_stages),
    has_isolation = .as_prob(has_isolation,    "has_isolation"),
    has_quarantine = .as_prob(has_quarantine,   "has_quarantine"),
    has_pep = .as_prob(has_pep,          "has_pep"),
    has_vaccination = .as_prob(has_vaccination,  "has_vaccination"),
    has_demography = .as_prob(has_demography,   "has_demography"),
    has_importation  = .as_prob(has_importation,  "has_importation")
  )

  # Checks
  # PEP acts on L and QL only; with no latent stage there is nothing to act on.
  if (s$has_pep > 0 && s$has_latent < 1) {
    if (s$has_latent == 0)
      stop("`has_pep = TRUE` requires `has_latent = TRUE` ",
           "(PEP acts on the L and QL compartments)")
    warning("`has_pep` is only active on runs where `has_latent` is drawn TRUE; ",
            "it will be forced off otherwise.")
  }

  # immunity = "none" makes the waning-stage count meaningless
  if (immunity == "none" && s$n_waning_stages != 1L) {
    warning('immunity = "none" (SIS) ignores `n_waning_stages`; forcing to 1')
    s$n_waning_stages <- 1L
  }

  # SI models have no recovery, so immunity fate is undefined
  if (s$has_recovery == 0 && immunity != "permanent") {
    warning("`has_recovery = FALSE` (SI model): immunity is undefined, ",
            'setting immunity = "permanent"')
    s$immunity <- "permanent"
  }

  structure(s, class = c("reservoir_structure", "list"))
}


#' @export
print.reservoir_structure <- function(x, ...) {
  fmt <- function(p) if (p >= 1) "TRUE" else if (p <= 0) "FALSE" else sprintf("p=%.2f", p)
  cat("<reservoir model structure>\n")
  cat(sprintf("  skeleton         : %s\n", .skeleton_name(x)))
  cat(sprintf("  has_latent       : %s   [STRUCTURAL]\n", fmt(x$has_latent)))
  cat(sprintf("  has_asymptomatic : %s\n", fmt(x$has_asymptomatic)))
  cat(sprintf("  has_recovery     : %s\n", fmt(x$has_recovery)))
  cat(sprintf("  immunity         : %s%s\n", x$immunity,
              if (x$immunity == "none") "   [STRUCTURAL]" else ""))
  cat(sprintf("  n_waning_stages  : %d\n", x$n_waning_stages))
  cat(sprintf("  has_isolation    : %s\n", fmt(x$has_isolation)))
  cat(sprintf("  has_quarantine   : %s\n", fmt(x$has_quarantine)))
  cat(sprintf("  has_pep          : %s\n", fmt(x$has_pep)))
  cat(sprintf("  has_vaccination  : %s\n", fmt(x$has_vaccination)))
  cat(sprintf("  has_demography   : %s\n", fmt(x$has_demography)))
  cat(sprintf("  has_importation  : %s\n", fmt(x$has_importation)))
  invisible(x)
}


################################################################################
# .skeleton_name - human-readable model name for a REALISED structure
# Uses L for the latent stage (package convention); SLIR == SEIR.
################################################################################
.skeleton_name <- function(s) {
  lat  <- s$has_latent       >= 1
  asy  <- s$has_asymptomatic >= 1
  rec  <- s$has_recovery     >= 1
  mixed <- any(vapply(s[c("has_latent","has_asymptomatic","has_recovery")],
                      function(p) p > 0 && p < 1, logical(1)))
  if (mixed) return("<mixture>")

  nm <- "S"
  if (lat) nm <- paste0(nm, "L")
  nm <- paste0(nm, "I")
  if (asy) nm <- paste0(nm, "A")
  if (!rec) return(nm) # SI / SLI / SLIA
  nm <- switch(s$immunity,
               none = paste0(nm, "S"),  # SIS  / SLIS
               permanent = paste0(nm, "R"),  # SIR  / SLIR
               waning  = paste0(nm, "RS")) # SIRS / SLIRS
  nm
}


################################################################################
# population_spec
#
# Population geometry: size, structure, horizon, mixing, seeding.
# NULL fields are sampled per run; supplied fields are pinned.
################################################################################
#' Specify population geometry
#'
#' @param pop_size Total population N. NULL samples log-uniformly from
#'   `pop_size_range`.
#' @param num_pops Number of patches P. NULL samples from `num_pops_probs`.
#' @param num_ages Number of age groups A. NULL samples from `num_ages_probs`.
#' @param num_days Simulation horizon. NULL samples from `num_days_probs`.
#' @param population_distribution Optional (P, A) matrix of population
#'   fractions. NULL generates a Dirichlet draw.
#' @param contact_matrix Optional (A, A) contact matrix. NULL generates a
#'   POLYMOD-style reciprocity-enforced matrix.
#' @param mixing_matrix Optional (P, P) patch-mixing matrix with rows summing
#'   to 1. NULL generates one.
#' @param seed_frac Fraction of N seeded as active infections. NULL draws
#'   log-uniformly from `seed_frac_range`.
#' @param novel TRUE for a fully susceptible population; FALSE seeds
#'   pre-existing immunity.
#' @param immune_frac Explicit immune fraction for non-novel runs. NULL derives
#'   it from the endemic equilibrium 1 - 1/R0.
#'
#' @export
population_spec <- function(
    pop_size = NULL,
    pop_size_range  = c(5e3, 4e7),

    num_pops = NULL,
    num_ages = NULL,
    num_pops_probs = c("1" = 0.67, "2" = 0.17, "3" = 0.16),
    num_ages_probs  = c("1" = 0.50, "2" = 0.17, "3" = 0.17, "5" = 0.16),

    num_days = NULL,
    num_days_probs = c("365" = 0.60, "730" = 0.20, "1095" = 0.15, "2000" = 0.05),

    population_distribution = NULL,
    contact_matrix = NULL,
    mixing_matrix = NULL,

    seed_frac  = NULL,
    seed_frac_range = c(1e-6, 5e-4),
    novel = TRUE,
    immune_frac = NULL
) {
  chk_probs <- function(p, nm) {
    if (any(p < 0)) stop("`", nm, "` must be non-negative")
    if (abs(sum(p) - 1) > 1e-8) stop("`", nm, "` must sum to 1")
    if (is.null(names(p))) stop("`", nm, "` must be a named vector, e.g. c('1'=0.5,'2'=0.5)")
    invisible(TRUE)
  }
  chk_probs(num_pops_probs, "num_pops_probs")
  chk_probs(num_ages_probs, "num_ages_probs")
  chk_probs(num_days_probs, "num_days_probs")

  if (!is.null(pop_size) && pop_size <= 0) stop("`pop_size` must be > 0")
  if (!is.null(num_days) && num_days <= 0) stop("`num_days` must be > 0")
  if (!is.null(num_pops) && num_pops <= 0) stop("`num_pops` must be > 0")
  if (!is.null(num_ages) && num_ages <= 0) stop("`num_ages` must be > 0")

  # Cross-checks between supplied matrices and supplied dimensions
  if (!is.null(contact_matrix)) {
    if (!is.matrix(contact_matrix) || nrow(contact_matrix) != ncol(contact_matrix))
      stop("`contact_matrix` must be a square A x A matrix")
    if (!is.null(num_ages) && nrow(contact_matrix) != num_ages)
      stop("`contact_matrix` is ", nrow(contact_matrix), "x", ncol(contact_matrix),
           " but `num_ages` = ", num_ages)
  }
  if (!is.null(mixing_matrix)) {
    if (!is.matrix(mixing_matrix) || nrow(mixing_matrix) != ncol(mixing_matrix))
      stop("`mixing_matrix` must be a square P x P matrix")
    if (!is.null(num_pops) && nrow(mixing_matrix) != num_pops)
      stop("`mixing_matrix` is ", nrow(mixing_matrix), "x", ncol(mixing_matrix),
           " but `num_pops` = ", num_pops)
    if (any(abs(rowSums(mixing_matrix) - 1) > 1e-8))
      stop("`mixing_matrix` rows must each sum to 1")
  }
  if (!is.null(population_distribution)) {
    pd <- as.matrix(population_distribution)
    if (!is.null(num_pops) && nrow(pd) != num_pops)
      stop("`population_distribution` has ", nrow(pd), " rows but `num_pops` = ", num_pops)
    if (!is.null(num_ages) && ncol(pd) != num_ages)
      stop("`population_distribution` has ", ncol(pd), " cols but `num_ages` = ", num_ages)
  }

  structure(list(
    pop_size = pop_size, pop_size_range = pop_size_range,
    num_pops = if (is.null(num_pops)) NULL else as.integer(num_pops),
    num_ages = if (is.null(num_ages)) NULL else as.integer(num_ages),
    num_pops_probs = num_pops_probs,
    num_ages_probs = num_ages_probs,
    num_days = if (is.null(num_days)) NULL else as.integer(num_days),
    num_days_probs = num_days_probs,
    population_distribution = population_distribution,
    contact_matrix = contact_matrix,
    mixing_matrix  = mixing_matrix,
    seed_frac = seed_frac, seed_frac_range = seed_frac_range,
    novel = novel, immune_frac = immune_frac
  ), class = c("reservoir_population", "list"))
}


#' @export
print.reservoir_population <- function(x, ...) {
  f <- function(v) if (is.null(v)) "<sampled>" else format(v)
  cat("<reservoir population spec>\n")
  cat(sprintf("  pop_size : %s\n", f(x$pop_size)))
  cat(sprintf("  num_pops : %s\n", f(x$num_pops)))
  cat(sprintf("  num_ages : %s\n", f(x$num_ages)))
  cat(sprintf("  num_days : %s\n", f(x$num_days)))
  cat(sprintf("  novel : %s\n", x$novel))
  invisible(x)
}


################################################################################
# intervention_policy
#
# Two modes:scheduled and reactive
################################################################################
#' Specify intervention policy
#'
#' @param mode "scheduled" (fixed calendar days) or "reactive"
#'   (threshold-triggered inside the simulation loop).
#' @param enable_prob Probability that a given run has any intervention at all.
#' @param start_day,duration Scheduled mode: NULL samples them.
#' @param on_threshold,off_threshold Reactive mode: daily case counts at which
#'   interventions engage and release. Given as a FRACTION of N if < 1,
#'   otherwise as absolute counts. NULL samples them.
#' @param trigger_delay Reactive mode: days between threshold crossing and
#'   interventions taking effect.
#' @param contact_reduction Reactive mode: multiplier applied to beta while
#'   interventions are active (school/work closure). 1.0 = no contact effect.
#' @param iso_rate,quar_rate,pep_rate,prep_start_rate,prep_eff Intervention
#'   intensities as daily rates (prep_eff is a 0-1 susceptibility reduction).
#'   NULL samples from the ranges in `reservoir_config()`.
#' @param school_closure,work_closure Contact reductions folded into beta.
#'
#' @export
intervention_policy <- function(
    mode          = c("scheduled", "reactive"),
    enable_prob   = 0.5,

    # scheduled
    start_day     = NULL,
    duration      = NULL,

    # reactive
    on_threshold  = NULL,
    off_threshold = NULL,
    trigger_delay = NULL,
    contact_reduction = NULL,

    # intensities (NULL = sample)
    iso_rate        = NULL,
    quar_rate       = NULL,
    pep_rate        = NULL,
    prep_start_rate = NULL,
    prep_eff        = NULL,

    # contact-reduction levers (folded into beta)
    school_closure  = NULL,
    work_closure    = NULL
) {
  mode <- match.arg(mode)
  enable_prob <- .as_prob(enable_prob, "enable_prob")

  if (!is.null(trigger_delay) && trigger_delay < 0)
    stop("`trigger_delay` must be >= 0")
  if (!is.null(on_threshold) && !is.null(off_threshold) &&
      off_threshold > on_threshold)
    stop("`off_threshold` must be <= `on_threshold` (otherwise the ",
         "controller oscillates every day)")

  structure(list(
    mode = mode, enable_prob = enable_prob,
    start_day = start_day, duration = duration,
    on_threshold = on_threshold, off_threshold = off_threshold,
    trigger_delay = trigger_delay, contact_reduction = contact_reduction,
    iso_rate = iso_rate, quar_rate = quar_rate, pep_rate = pep_rate,
    prep_start_rate = prep_start_rate, prep_eff = prep_eff,
    school_closure = school_closure, work_closure = work_closure
  ), class = c("reservoir_policy", "list"))
}


#' @export
print.reservoir_policy <- function(x, ...) {
  cat("<reservoir intervention policy>\n")
  cat(sprintf("  mode        : %s\n", x$mode))
  cat(sprintf("  enable_prob : %.2f\n", x$enable_prob))
  invisible(x)
}


################################################################################
# reservoir_config
#
# Priors, observation model, and output scope.
# sample_params.R.
################################################################################
#' Configure priors, observation model, and output scope
#'
#' @export
reservoir_config <- function(
    # R0 priors, per mechanism
    r0_base = list(meanlog = log(3.0), sdlog = 0.7, shift = 1.0, floor = 1.1),
    r0_vector = list(meanlog = log(2.5), sdlog = 0.7, shift = 0.5, floor = 1.0),
    r0_waterborne = list(meanlog = log(2.5), sdlog = 0.6, shift = 0.5, floor = 1.0),

    # epidemiological params
    infectious_days = list(meanlog = log(6), sdlog = 0.8, lo = 0.5, hi = 40),
    latent_days = list(meanlog = log(6), sdlog = 0.8, lo = 0.5, hi = 30),
    immunity_days = list(meanlog = log(365), sdlog = 1.2, lo = 7, hi = 3650 * 5),
    p_asym_beta = c(1.5, 4),
    asymp_trans_beta = c(2, 5),
    kappa_logrange = c(log(0.02), log(0.5)),
    beta_bounds = c(1e-4, 5.0),

    # direct beta (required when has_recovery = FALSE)
    beta_direct = NULL,

    # waves
    wave_count_probs = c("0" = 0.35, "1" = 0.30, "2" = 0.20, "3" = 0.10, "4" = 0.05),
    wave_mult_sdlog = 0.6,
    wave_crash_prob = 0.30, wave_crash_range = c(0.05, 0.25),
    wave_surge_prob = 0.20, wave_surge_range = c(1.5, 3.0),
    wave_dipsurge_prob = 0.25,

    # seasonality
    seasonality_prob = 0.8,
    n_harmonics_range = c(1L, 3L),
    seasonality_amp_range = c(0.1, 0.5),
    annual_jitter_max = 30,
    daily_noise_sd = 0.05,
    daily_noise_clamp  = c(0.5, 2.0),

    # super-spreading
    # Super-spreading.
    superspread_event_prob_range = c(0.005, 0.10),   # daily P(event)
    superspread_boost_range = c(0.5, 4.0),      # mean excess on event days
    superspread_k_logrange = c(log(0.5), log(3.0)),  # dispersion of excess

    # demography / importation
    birth_rate_range = c(2e-5, 1.2e-4),
    death_rate_factor = c(0.8, 1.3),
    importation_logrange = c(log(0.01), log(5.0)),

    # intervention intensity ranges
    iso_rate_range = c(0.05, 0.5),
    quar_rate_range = c(0.02, 0.3),
    pep_rate_range = c(0.01, 0.2),
    prep_start_rate_range = c(0.001, 0.02),
    prep_eff_range = c(0.5, 0.95),
    school_closure_range  = c(0.3, 0.8),
    work_closure_range = c(0.3, 0.8),
    interv_start_frac = c(0.02, 0.5),
    interv_duration_range = c(14L, 180L),

    # reactive controller ranges
    on_threshold_frac_range = c(5e-5, 1e-3),
    off_threshold_frac_range = c(0.0, 5e-4),
    trigger_delay_range = c(0L, 14L),
    contact_reduction_range  = c(0.3, 0.9),

    # clinical severity
    hosp_prob_logrange = c(log(0.001), log(0.90)),
    death_prob_logrange = c(log(0.001), log(0.75)),

    # observation model
    detection_logrange = c(log(0.01), log(0.60)),
    report_delay_range = c(2, 21),
    r_nb = list(meanlog = log(5), sdlog = 0.8),
    underreport_prob = 0.8,
    weekday_prob = 0.8,
    lab_batch_prob = 0.8,
    report_delay_prob = 0.8,

    # mechanism-specific
    vector = list(
      sigma_v_logrange = c(log(1/21), log(1/4)),
      mu_v_logrange = c(log(1/30), log(1/7)),
      a_bite_range = c(0.1, 1.0),
      b_h_range = c(0.1, 0.8),
      b_v_range = c(0.1, 0.8),
      Nv_init_frac_range = c(0.05, 2.0),
      Iv_frac_range = c(0.001, 0.05)
    ),
    waterborne = list(
      eta_I_logrange = c(log(0.1), log(5.0)),
      eta_A_logrange = c(log(0.05), log(2.0)),
      alpha_env_range = c(0.1, 1.0),
      mu_w_logrange = c(log(0.01), log(0.5)),
      waterborne_frac_range = c(0.1, 0.9),
      delta_eff_bounds = c(1e-14, 1.0)
    ),

    # tau-leaping controls stuff
    epsilon = 0.03, Ncritical = 10L, exactThreshold = 10.0, maxtau = Inf,

    # convenience
    # TRUE switches off every source of time-varying / stochastic forcing
    # (seasonality, waves, daily noise, super-spreading).
    deterministic       = FALSE,

    # output signals
    compute_rt          = TRUE,
    compute_noise_full  = TRUE,
    compute_wastewater  = TRUE,
    compute_syndromic   = TRUE,
    compute_hosp_deaths = TRUE
) {
  cfg <- as.list(environment())

  # Distribution fields accept a scalar, a length-2 range, or a full list.
  for (nm in c("r0_base", "r0_vector", "r0_waterborne"))
    cfg[[nm]] <- .as_dist(cfg[[nm]], nm, "shifted")
  for (nm in c("infectious_days", "latent_days", "immunity_days", "r_nb"))
    cfg[[nm]] <- .as_dist(cfg[[nm]], nm, "clipped")

  if (isTRUE(deterministic)) {
    cfg$seasonality_prob             <- 0
    cfg$wave_count_probs             <- c("0" = 1)
    cfg$daily_noise_sd               <- 0
    cfg$superspread_event_prob_range <- c(0, 0)
    cfg$superspread_boost_range      <- c(0, 0)
  }
  cfg$deterministic <- NULL

  structure(cfg, class = c("reservoir_config", "list"))
}


#' @export
print.reservoir_config <- function(x, ...) {
  cat("<reservoir config>\n")
  cat(sprintf("  R0 (base)     : LogNormal(%.2f, %.2f) + %.1f\n",
              x$r0_base$meanlog, x$r0_base$sdlog, x$r0_base$shift))
  cat(sprintf("  seasonality   : p = %.2f\n", x$seasonality_prob))
  cat(sprintf("  output scope  : rt=%s noise_full=%s ww=%s synd=%s hosp=%s\n",
              x$compute_rt, x$compute_noise_full, x$compute_wastewater,
              x$compute_syndromic, x$compute_hosp_deaths))
  invisible(x)
}


################################################################################
# .scope_config - apply a coarse `scope` argument on top of a config
#
#   "full"       everything
#   "signals"    all observation streams but no Rt (skips the NGM loop)
#   "cases_only" true + limited-noise case counts only (fastest)
################################################################################
.scope_config <- function(config, scope) {
  scope <- match.arg(scope, c("full", "signals", "cases_only"))
  if (scope == "full") return(config)
  if (scope == "signals") {
    config$compute_rt <- FALSE
    return(config)
  }
  config$compute_rt  <- FALSE
  config$compute_noise_full <- FALSE
  config$compute_wastewater <- FALSE
  config$compute_syndromic <- FALSE
  config$compute_hosp_deaths <- FALSE
  config
}


################################################################################
# .sample_structure - draw one concrete structure from a (possibly
# probabilistic) model_structure(). Returns a list of hard TRUE/FALSE.
################################################################################
.sample_structure <- function(s) {
  out <- list(
    has_latent       = .draw_flag(s$has_latent),
    has_asymptomatic = .draw_flag(s$has_asymptomatic),
    has_recovery     = .draw_flag(s$has_recovery),
    immunity         = s$immunity,
    n_waning_stages  = s$n_waning_stages,
    has_isolation    = .draw_flag(s$has_isolation),
    has_quarantine   = .draw_flag(s$has_quarantine),
    has_pep          = .draw_flag(s$has_pep),
    has_vaccination  = .draw_flag(s$has_vaccination),
    has_demography   = .draw_flag(s$has_demography),
    has_importation  = .draw_flag(s$has_importation)
  )
  # PEP needs a latent stage to act on
  if (!out$has_latent) out$has_pep <- FALSE
  # SIS ignores waning stages
  if (out$immunity == "none") out$n_waning_stages <- 1L
  out$skeleton <- .skeleton_name(out)
  out
}
