################################################################################
# generate_dataset.R
#
#  draw a structure, draw parameters, build a cfg, run the
#  simulators, post-process, label, and collect.
#
################################################################################


################################################################################
# cfg builders
################################################################################

# Build a configuration object for human-to-human
#' @export
build_cfg_base <- function(p) {

  s <- p$structure

  cfg <- setup_sim(
    N = p$N, P = p$P, A = p$A, C = p$C, Mp = p$Mp,
    pop_fracs = p$pop_fracs,
    tf_days = p$tf_days,
    beta = p$beta,
    infectious_days = p$infectious_days,
    latent_days = p$latent_days,
    immunity_days = p$immunity_days,
    p_asym = p$p_asym,
    asymp_trans = p$asymp_trans,
    kappa = p$kappa,
    n_waning_stages = p$n_waning_stages,

    # structural flags
    no_latent = !s$has_latent,
    sis_mode = (s$immunity == "none"),

    # seeding is overwritten by set_initial_conditions(); seed into a
    # compartment that is guaranteed to exist under the chosen structure
    n_seed = 1,
    seed_comp = "Is",

    #Demographics params
    birth_rate = p$birth_rate,
    death_rate = p$death_rate,
    importation_rate = p$importation_rate,

    #Wave change params
    wave_change_days = p$wave_change_days,
    wave_change_betas = p$wave_change_betas,

    #Seasonality and noisiness params
    use_seasonality  = p$use_seasonality,
    n_harmonics = p$n_harmonics,
    harmonic_amps = p$harmonic_amps,
    harmonic_offsets = p$harmonic_offsets,
    harmonic_periods = p$harmonic_periods,
    annual_jitter_sd = p$annual_jitter_sd,
    daily_noise_sd = p$daily_noise_sd,
    daily_noise_clamp = p$daily_noise_clamp %||% c(0.5, 2.0),

    # school / work closure fold into beta
    contact_mult = (1 - p$school_closure) * (1 - p$work_closure),

    iso_on_day = p$iso_on_day,  iso_off_day  = p$iso_off_day,  iso_active  = p$iso_active,
    quar_on_day = p$quar_on_day, quar_off_day = p$quar_off_day, quar_active = p$quar_active,
    pep_on_day = p$pep_on_day,  pep_off_day  = p$pep_off_day,  pep_active  = p$pep_active,
    prep_on_day = p$prep_on_day, prep_off_day = p$prep_off_day, prep_active = p$prep_active,

    # PrEP/vaccination campaign rate. Without this, Sprep never populates
    # and prep_active has no effect on the trajectory.
    prep_start_rate = p$prep_start_rate,

    reactive = identical(p$interv_mode, "reactive") && p$has_intervention,
    on_threshold = p$on_threshold,
    off_threshold = p$off_threshold,
    trigger_delay = p$trigger_delay,
    reactive_contact_mult = p$contact_reduction
  )

  # Initial conditions
  cfg <- set_initial_conditions(
    cfg,
    novel = p$novel,
    immune_frac = p$immune_frac,
    R0 = p$R0,
    seed_frac = p$seed_frac,
    infectious_days  = if (is.finite(p$infectious_days)) p$infectious_days else 10,
    latent_days = p$latent_days,
    p_asym = p$p_asym,
    has_intervention = p$has_intervention,
    waning_stages = p$n_waning_stages
  )

  # Super-spreading
  # Daily stochastic multiplier on beta; normalized so E[multiplier] = 1 to
  # preserve mean transmission.
  if (!is.null(p$ss_event_prob) && p$ss_event_prob > 0 && p$ss_boost > 0) {
    tf <- cfg$tf_days
    event_prob <- pmin(1, p$ss_event_prob)

    # Excess is Gamma(k, boost/k): mean = boost, always >= 0.
    is_event <- runif(tf) < event_prob
    excess <- ifelse(is_event,
                       rgamma(tf, shape = p$ss_shape, scale = p$ss_scale), 0.0)
    multipliers <- 1 + excess

    # Rescale so E[multiplier] = 1, preserving expected transmission.
    expected <- 1 + event_prob * p$ss_boost
    multipliers <- multipliers / expected

    cfg$beta_eff_daily <- pmax(0, cfg$beta_eff_daily * multipliers)
  }

  # Stash structure
  cfg$structure <- s
  cfg
}


#' Build a configuration object for the vector-borne engine
#' @export
build_cfg_vector <- function(p) {
  cfg <- build_cfg_base(p)

  Nv <- round(p$N * p$Nv_init_frac)
  n_Iv <- round(Nv * p$Iv_frac)
  n_Sv <- Nv - n_Iv
  cfg$x0 <- c(cfg$x0, n_Sv, 0L, n_Iv, 0L)

  cfg$params$sigma_v <- p$sigma_v
  cfg$params$mu_v <- p$mu_v
  cfg$params$a_bite <- p$a_bite
  cfg$params$b_h <- p$b_h
  cfg$params$b_v <- p$b_v

  cfg$seasonal_daily <- rep(1.0, p$tf_days)
  cfg$f_interv_daily <- rep(1.0, p$tf_days)

  cfg$comp_names <- c(.HUMAN_COMPS, "Sv", "Ev", "Iv", "cumul_v")
  cfg
}


#' Build a configuration object for the environmental engine
#' @export
build_cfg_waterborne <- function(p) {
  cfg <- build_cfg_base(p)

  cfg$x0 <- c(cfg$x0, 0.0)

  cfg$params$eta_I <- p$eta_I
  cfg$params$eta_A <- p$eta_A
  cfg$params$alpha_env <- p$alpha_env
  cfg$params$mu_w <- p$mu_w

  cfg$delta_eff_daily <- rep(p$delta_eff, p$tf_days)
  cfg$f_contact_daily <- rep(1.0, p$tf_days)
  cfg$f_water_daily <- rep(1.0, p$tf_days)
  cfg$seasonal_daily <- rep(1.0, p$tf_days)

  cfg$comp_names <- c(.HUMAN_COMPS, "W")
  cfg
}


################################################################################
# run_one: one simulation, from parameters to labeled output
################################################################################
run_one <- function(mechanism, params, run_id = 1L,
                    config = reservoir_config(),
                    min_attack_rate = 0.0) {

  cfg <- tryCatch(
    switch(mechanism,
           base = build_cfg_base(params),
           vector = build_cfg_vector(params),
           waterborne = build_cfg_waterborne(params)),
    error = function(e) {
      message("Run", run_id, "cfg error:", conditionMessage(e)); NULL
    })
  if (is.null(cfg)) return(NULL)

  result <- tryCatch(.dispatch_engine(mechanism, cfg),
    error = function(e) {
      message("Run", run_id, "failed:", conditionMessage(e)); NULL
    })
  if (is.null(result)) return(NULL)

  daily_true <- rowSums(result$daily_new_cases_block)

  # Die-out filter
  total <- sum(daily_true)
  if (min_attack_rate > 0) {
    if (total / max(cfg$N, 1) < min_attack_rate) {
      message("Run", run_id, "-attack rate",
              signif(total / max(cfg$N, 1), 3), "below", min_attack_rate,
              ", skipping")
      return(NULL)
    }
  } else if (total < 1) {
    message("Run", run_id, "- no cases, skipping")
    return(NULL)
  }

  # Limited noise (always computed; cheap)
  nl <- apply_noise_limited(daily_true)

  # Full observation pipeline (optional)
  nf <- NULL
  if (isTRUE(config$compute_noise_full)) {
    nf <- tryCatch(
      run_postprocessing(
        result = result, params = cfg$params,
        P = cfg$P, A = cfg$A, num_days = cfg$tf_days,
        age_hosp_probs  = params$age_hosp_probs,
        age_death_probs = params$age_death_probs,
        config = config),
      error = function(e) {
        message("Run", run_id, "postprocessing error:", conditionMessage(e))
        NULL
      })
  }

  # Rt (optional; the NGM loop is expensive)
  Rt <- if (isTRUE(config$compute_rt)) {
    tryCatch(compute_daily_rt(result, cfg, mechanism = mechanism),
             error = function(e) rep(NA_real_, cfg$tf_days))
  } else rep(NA_real_, cfg$tf_days)

  nz_Rt <- Rt[is.finite(Rt) & Rt > 0]
  s <- params$structure

  summ <- list(
    total_cases = total,
    infections_per_person = round(total / max(cfg$N, 1), 4),
    attack_rate = round(total / max(cfg$N, 1), 4),
    peak_day = which.max(daily_true),
    peak_cases = max(daily_true),
    peak_Rt = if (length(nz_Rt)) round(max(nz_Rt), 3) else NA,
    final_Rt = if (length(nz_Rt)) round(tail(nz_Rt, 1), 3) else NA,
    total_hosp = if (!is.null(nf)) sum(nf$daily_hosp_true)  else NA,
    total_deaths = if (!is.null(nf)) sum(nf$daily_death_true) else NA,
    R0  = params$R0,
    mechanism = mechanism,
    # structure modular labels
    skeleton = s$skeleton,
    has_latent = s$has_latent,
    has_asymptomatic  = s$has_asymptomatic,
    has_recovery = s$has_recovery,
    immunity = s$immunity,
    n_waning_stages = s$n_waning_stages,
    has_isolation = s$has_isolation,
    has_quarantine = s$has_quarantine,
    has_pep = s$has_pep,
    has_vaccination = s$has_vaccination,
    has_demography = s$has_demography,
    has_importation = s$has_importation,
    has_intervention = params$has_intervention,
    interv_mode = params$interv_mode,
    is_vectorborne = as.integer(mechanism == "vector"),
    is_waterborne = as.integer(mechanism == "waterborne")
  )

  list(
    run_id = run_id,
    mechanism = mechanism,
    scenario_type = if (params$endemic) "endemic" else "nonendemic",
    params = params,
    cfg = cfg,
    daily_true = daily_true,
    daily_Rt = Rt,
    noise_limited = nl,
    noise_full = nf,
    state_over_time = result$state_over_time,
    intervention_active = result$intervention_active,
    final_state = result$final_state,
    cumul = result$cumul_infections,
    summary = summ
  )
}


################################################################################
# generate_dataset
################################################################################
#' Generate a labelled synthetic epidemic dataset
#'
#' @param n_per_mechanism Number of valid runs to collect per mechanism.
#' @param mechanisms Any of "base", "vector", "waterborne".
#' @param structure A `model_structure()`. Fields may be probabilities, in
#'   which case the structure is drawn per run and labelled in the output.
#' @param population A `population_spec()`.
#' @param policy An `intervention_policy()`.
#' @param config A `reservoir_config()`.
#' @param scope "full", "signals" (no Rt), or "cases_only" (fastest).
#'   Overrides the individual compute_* flags in `config`.
#' @param min_attack_rate Discard runs below this cumulative attack rate.
#' @param max_attempts Give up after this many draws per mechanism.
#' @param save_dir Directory for per-mechanism .rds files. NULL keeps results
#'   in memory only.
#' @param seed Base RNG seed; run i uses `seed + i` for reproducibility.
#'
#' @export
generate_dataset <- function(
    n_per_mechanism = 50,
    mechanisms      = c("base", "vector", "waterborne"),
    structure       = model_structure(),
    population      = population_spec(),
    policy          = intervention_policy(),
    config          = reservoir_config(),
    scope           = c("full", "signals", "cases_only"),
    min_attack_rate = 0.0,
    max_attempts    = 200,
    save_dir        = NULL,
    seed            = NULL,
    verbose         = TRUE
) {

  scope  <- match.arg(scope)
  config <- .scope_config(config, scope)

  if (!inherits(structure, "reservoir_structure"))
    stop("`structure` must come from model_structure()")
  if (!inherits(population, "reservoir_population"))
    stop("`population` must come from population_spec()")
  if (!inherits(policy, "reservoir_policy"))
    stop("`policy` must come from intervention_policy()")

  dataset <- list()

  for (mech in mechanisms) {

    if (verbose) message("\n=== Mechanism: ", toupper(mech),
                         "  (scope: ", scope, ") ===")

    runs <- list(); attempts <- 0L; run_id <- 1L

    while (length(runs) < n_per_mechanism && attempts < max_attempts) {
      attempts <- attempts + 1L
      if (!is.null(seed)) set.seed(seed + attempts)

      if (verbose) message("  [", mech, "] attempt ", attempts,
                           " | collected ", length(runs), "/", n_per_mechanism)

      params <- tryCatch(
        sample_params(mech, structure, population, policy, config),
        error = function(e) {
          message("  Sampler failed: ", conditionMessage(e)); NULL
        })
      if (is.null(params)) next

      run <- run_one(mech, params, run_id = run_id,
                     config = config, min_attack_rate = min_attack_rate)

      if (!is.null(run)) {
        runs[[length(runs) + 1L]] <- run
        run_id <- run_id + 1L
      }
    }

    if (length(runs) < n_per_mechanism)
      warning("Only collected ", length(runs), "/", n_per_mechanism,
              " runs for '", mech, "' after ", max_attempts, " attempts")

    dataset[[mech]] <- runs

    if (!is.null(save_dir)) {
      dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
      out_path <- file.path(save_dir, paste0(mech, "_runs.rds"))
      saveRDS(runs, out_path)
      if (verbose) message("  Saved to: ", out_path)
    }
  }

  dataset
}

