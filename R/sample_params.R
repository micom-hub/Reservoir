################################################################################
# sample_params.R
#
# hierarchical sampler
#
################################################################################

################################################################################
# small helpers
################################################################################
.runif_r  <- function(r) runif(1, r[1], r[2])
.logunif  <- function(r) exp(runif(1, r[1], r[2]))
.rlnorm_c <- function(p) pmax(p$lo, pmin(rlnorm(1, p$meanlog, p$sdlog), p$hi))

# Draw an R0 from a distribrition
.draw_r0 <- function(p) {
  raw <- if (!is.null(p$shape))
           rgamma(1, shape = p$shape, rate = p$rate %||% 1) #gamma
         else
           rlnorm(1, p$meanlog, p$sdlog) #lognorm dist
  max(raw + (p$shift %||% 0), p$floor %||% 0)
}

.sample_named <- function(probs) {
  .sample1(as.numeric(names(probs)), prob = probs)
}


################################################################################
# compute_ip (infectious poetnital)
################################################################################
compute_ip <- function(infectious_days, latent_days, p_asym,
                       asymp_trans, kappa, has_latent = TRUE) {
  gamma <- 1 / max(infectious_days, 1e-8)
  sigma <- 1 / max(latent_days,     1e-8)

  mu_Is <- gamma #Is and Ia are symp and asymp
  mu_Ia <- kappa + gamma

  tau_Is <- 1 / max(mu_Is, 1e-8)
  tau_Ia <- 1 / max(mu_Ia, 1e-8)

  # Asymptomatic contribute
  contrib_Ia <- asymp_trans * tau_Ia + (kappa / max(mu_Ia, 1e-8)) * tau_Is

  if (has_latent) {
    mu_L <- sigma
    ip <- (sigma * (1 - p_asym) / max(mu_L, 1e-8)) * tau_Is +
          (sigma * p_asym       / max(mu_L, 1e-8)) * contrib_Ia
  } else {
    # No latent stage 
    ip <- (1 - p_asym) * tau_Is + p_asym * contrib_Ia
  }

  max(ip, 1e-8)
}


################################################################################
# sample_structural_multiplier
#
################################################################################
sample_structural_multiplier <- function(P, A, C, Mp, pop_fracs) {
  if (P == 1 && A == 1) return(1.0)
  N_flat <- rep(rep(1/P, P), each = A) * rep(pop_fracs, times = P)
  inv_N  <- ifelse(N_flat > 0, 1 / N_flat, 0)
  M_raw  <- kronecker(Mp, C)
  M_0   <- M_raw * outer(inv_N, N_flat)
  max(Re(eigen(M_0, only.values = TRUE)$values))
}


################################################################################
# sample_contact_matrix- POLYMOD-like symmetric contact matrix
################################################################################
sample_contact_matrix <- function(A, pop_fracs = NULL) {
  if (A == 1) return(matrix(1.0, 1, 1))
  if (is.null(pop_fracs)) pop_fracs <- rep(1/A, A)

  C_raw <- matrix(0, A, A)
  for (i in seq_len(A))
    for (j in seq_len(A))
      C_raw[i, j] <- exp(-abs(i - j) * runif(1, 0.3, 0.8))

  if (A >= 2) {
    boost <- runif(1, 0.1, 0.3)
    C_raw[1, 2] <- C_raw[1, 2] + boost
    C_raw[2, 1] <- C_raw[2, 1] + boost
  }
  diag(C_raw) <- diag(C_raw) * runif(A, 1.5, 3.0)

  N <- pop_fracs
  C_sym <- matrix(0, A, A)
  for (i in seq_len(A))
    for (j in seq_len(A))
      C_sym[i, j] <- (N[i] * C_raw[i, j] + N[j] * C_raw[j, i]) /
                     (2 * max(N[i], 1e-8))
  C_sym / max(C_sym)
}


################################################################################
# sample_mixing_matrix — patch mixing (P x P), rows sum to 1
################################################################################
sample_mixing_matrix <- function(P) {
  if (P == 1) return(matrix(1.0, 1, 1))
  coupling <- runif(1, 0.01, 0.30)
  Mp <- matrix(coupling / (P - 1), P, P)
  diag(Mp) <- 1 - coupling
  Mp / rowSums(Mp)
}


################################################################################
# sample_population_structure
#
# Matches default values in a population_spec(); samples the rest.
################################################################################
sample_population_structure <- function(population) {
  A <- if (!is.null(population$num_ages)) population$num_ages
       else as.integer(.sample_named(population$num_ages_probs))
  P <- if (!is.null(population$num_pops)) population$num_pops
       else as.integer(.sample_named(population$num_pops_probs))

  # Age fractions
  if (!is.null(population$population_distribution)) {
    pd <- as.matrix(population$population_distribution)
    if (nrow(pd) != P || ncol(pd) != A)
      stop("`population_distribution` is ", nrow(pd), "x", ncol(pd),
           " but sampled P=", P, ", A=", A,
           ". Pin `num_pops`/`num_ages` when supplying a distribution.")
    pop_fracs <- colSums(pd) / sum(pd)
  } else if (A == 1) {
    pop_fracs <- 1.0
  } else {
    alpha<- runif(1, 0.5, 5.0)
    raw <- rgamma(A, shape = alpha, rate = 1)
    pop_fracs <- raw / sum(raw)
  }

  # Age multipliers for hospitalzation / death that increase w/age
  age_multipliers <- exp(seq(log(0.2), log(3.0), length.out = A))
  age_multipliers <- age_multipliers / mean(age_multipliers)

  C <- if (!is.null(population$contact_matrix)) {
         cm <- as.matrix(population$contact_matrix)
         if (nrow(cm) != A)
           stop("`contact_matrix` is ", nrow(cm), "x", ncol(cm),
                " but sampled A=", A, ". Pin `num_ages` when supplying one.")
         cm / max(cm)
       } else sample_contact_matrix(A, pop_fracs)

  Mp <- if (!is.null(population$mixing_matrix)) {
          mm <- as.matrix(population$mixing_matrix)
          if (nrow(mm) != P)
            stop("`mixing_matrix` is ", nrow(mm), "x", ncol(mm),
                 " but sampled P=", P, ". Pin `num_pops` when supplying one.")
          mm
        } else sample_mixing_matrix(P)

  list(P = P, A = A, C = C, Mp = Mp,
       pop_fracs = pop_fracs, age_multipliers = age_multipliers,
       structural_mult = sample_structural_multiplier(P, A, C, Mp, pop_fracs))
}


################################################################################
# sample_shared_epi
#

################################################################################
sample_shared_epi <- function(structure, config) {

  infectious_days <- .rlnorm_c(config$infectious_days)
  latent_days     <- .rlnorm_c(config$latent_days)
  immunity_days   <- .rlnorm_c(config$immunity_days)

  if (structure$has_asymptomatic) {
    p_asym      <- rbeta(1, config$p_asym_beta[1], config$p_asym_beta[2])
    asymp_trans <- rbeta(1, config$asymp_trans_beta[1], config$asymp_trans_beta[2])
    kappa       <- .logunif(config$kappa_logrange)
  } else {
    p_asym      <- 0.0
    asymp_trans <- 0.0
    kappa       <- 0.0
  }

  # Immunity definitions. "none" (SIS) is handled structurally via sis_mode; omega
  # is set to 0 so the still-present R -> S reaction stays inert.
  omega_days <- if (structure$immunity == "waning") immunity_days else Inf

  # Recovery. gamma = 0 gives SI / SLI.
  if (!structure$has_recovery) infectious_days <- Inf

  list(
    infectious_days = infectious_days,
    latent_days     = latent_days,
    immunity_days   = omega_days,
    p_asym  = p_asym,
    asymp_trans = asymp_trans,
    kappa  = kappa,
    n_waning_stages = structure$n_waning_stages
  )
}


################################################################################
# sample_scenario
#
# Waves, seasonality, super-spreading, interventions, demographics,
# importation, seeding, severity, observation model.
################################################################################
sample_scenario <- function(beta, structure, population, policy, config) {

  tf_days <- if (!is.null(population$num_days)) population$num_days
             else as.integer(.sample_named(population$num_days_probs))

  N <- if (!is.null(population$pop_size)) as.numeric(population$pop_size)
       else round(exp(runif(1, log(population$pop_size_range[1]),
                               log(population$pop_size_range[2]))))

  # Waves
  n_waves <- as.integer(.sample_named(config$wave_count_probs))
  wave_change_days <- NULL; wave_change_betas <- NULL

  if (n_waves > 0 && tf_days > 60) {
    # The pool can be smaller than n_waves on short horizons
    pool <- seq(50L, max(51L, tf_days - 10L))
    if (n_waves > length(pool)) n_waves <- length(pool)
    wave_change_days <- sort(if (n_waves == 1L) .sample1(pool)
                             else sample(pool, n_waves))
    wave_mult <- exp(rnorm(n_waves, 0, config$wave_mult_sdlog))

    if (runif(1) < config$wave_crash_prob)
      wave_mult[sample(n_waves, 1)] <- .runif_r(config$wave_crash_range)
    if (runif(1) < config$wave_surge_prob)
      wave_mult[sample(n_waves, 1)] <- .runif_r(config$wave_surge_range)
    if (n_waves >= 2 && runif(1) < config$wave_dipsurge_prob) {
      k <- .sample1(seq_len(n_waves - 1))
      wave_mult[k]     <- runif(1, 0.3, 0.6)
      wave_mult[k + 1] <- runif(1, 1.4, 2.5)
    }
    wave_change_betas <- pmax(config$beta_bounds[1],
                              pmin(beta * wave_mult, config$beta_bounds[2]))
  }

  # Seasonality
  use_seasonality <- runif(1) < config$seasonality_prob
  if (use_seasonality) {
    n_harmonics <- .sample1(config$n_harmonics_range[1]:config$n_harmonics_range[2])
    base_amp  <- .runif_r(config$seasonality_amp_range)
    fracs <- runif(n_harmonics, 0.3, 1.0); fracs <- fracs / sum(fracs)
    harmonic_amps <- base_amp * fracs
    harmonic_offsets <- runif(n_harmonics, 0, 365)
    harmonic_periods <- c(365.0, 182.5, 91.25)[seq_len(min(n_harmonics, 3))]
    annual_jitter_sd <- runif(1, 0, config$annual_jitter_max)
  } else {
    n_harmonics <- 0L; harmonic_amps <- numeric(0)
    harmonic_offsets <- numeric(0); harmonic_periods <- numeric(0)
    annual_jitter_sd <- 0
  }

  # Day-to-day environmental stochasticity on beta
  daily_noise_sd <- config$daily_noise_sd

  # Super-spreading
  # ss_event_prob : daily probability of a super-spreading event
  # ss_boost : mean multiplicative excess on an event day
  # ss_shape (k)  : dispersion of that excess (lower = burstier)
  ss_event_prob <- .runif_r(config$superspread_event_prob_range)
  ss_boost <- .runif_r(config$superspread_boost_range)
  ss_shape  <- .logunif(config$superspread_k_logrange)
  ss_scale  <- ss_boost / max(ss_shape, 1e-3)

  # Interventions
  any_interv <- (structure$has_isolation || structure$has_quarantine ||
                 structure$has_pep || structure$has_vaccination)
  has_intervention <- any_interv && (runif(1) < policy$enable_prob)

  pick <- function(fixed, rng) if (!is.null(fixed)) fixed else .runif_r(rng)

  iso_active <- 0; quar_active <- 0; pep_active <- 0
  prep_start <- 0; prep_eff <- 0
  school_red <- 0; work_red <- 0
  start_day <- NULL; end_day <- NULL
  on_threshold <- 0; off_threshold <- 0; trigger_delay <- 0L
  contact_reduction <- 1.0

  if (has_intervention) {
    if (structure$has_isolation)
      iso_active <- pick(policy$iso_rate, config$iso_rate_range)
    if (structure$has_quarantine)
      quar_active <- pick(policy$quar_rate, config$quar_rate_range)
    if (structure$has_pep)
      pep_active <- pick(policy$pep_rate, config$pep_rate_range)
    if (structure$has_vaccination) {
      prep_start <- pick(policy$prep_start_rate, config$prep_start_rate_range)
      prep_eff   <- pick(policy$prep_eff,        config$prep_eff_range)
    }

    if (policy$mode == "scheduled") {
      # Contact reduction is a scheduled measure: in reactive mode the engine
      # applies contact_reduction itself, so these stay at 0 there.
      school_red <- policy$school_closure %||% .runif_r(config$school_closure_range)
      work_red   <- policy$work_closure   %||% .runif_r(config$work_closure_range)

      lo <- max(10L, round(tf_days * config$interv_start_frac[1]))
      hi <- max(lo + 1L, round(tf_days * config$interv_start_frac[2]))
      start_day <- policy$start_day %||% .sample1(seq(lo, hi))
      dur <- policy$duration %||%
             .sample1(config$interv_duration_range[1]:config$interv_duration_range[2])
      end_day <- min(start_day + dur, tf_days)
    } else {
      # reactive: thresholds given as fractions of N when < 1
      ont <- policy$on_threshold  %||% .runif_r(config$on_threshold_frac_range)
      oft <- policy$off_threshold %||% .runif_r(config$off_threshold_frac_range)
      on_threshold  <- if (ont < 1) ont * N else ont
      off_threshold <- if (oft < 1) oft * N else oft
      if (off_threshold > on_threshold) off_threshold <- on_threshold
      trigger_delay <- as.integer(policy$trigger_delay %||%
        .sample1(config$trigger_delay_range[1]:config$trigger_delay_range[2]))
      contact_reduction <- policy$contact_reduction %||%
        .runif_r(config$contact_reduction_range)
    }
  }

  # Demographics
  if (structure$has_demography) {
    birth_rate <- .runif_r(config$birth_rate_range)
    death_rate <- birth_rate * .runif_r(config$death_rate_factor)
  } else {
    birth_rate <- 0.0; death_rate <- 0.0
  }

  importation_rate <- if (structure$has_importation)
    .logunif(config$importation_logrange) else 0.0

  seed_frac <- population$seed_frac %||%
    exp(runif(1, log(population$seed_frac_range[1]),
                 log(population$seed_frac_range[2])))

  # Count-mode seeding: an absolute number of seeds placed in one compartment
  seed_mode <- population$seed_mode %||% "fraction"
  n_seed <- if (identical(seed_mode, "count"))
    .sample1(seq.int(population$n_seed_range[1], population$n_seed_range[2]))
    else NA_integer_
  seed_comp <- if (identical(seed_mode, "count"))
    .sample1(names(population$seed_comp_probs),
             prob = population$seed_comp_probs)
    else NA_character_

  hosp_prob  <- .runif_r(config$hosp_prob_range)
  death_prob <- .runif_r(config$death_prob_range)
  detection <- .logunif(config$detection_logrange)
  report_delay <- round(.runif_r(config$report_delay_range))
  r_nb <- rlnorm(1, config$r_nb$meanlog, config$r_nb$sdlog)

  list(
    N = N, tf_days = tf_days,
    wave_change_days = wave_change_days, wave_change_betas = wave_change_betas,
    use_seasonality = use_seasonality, n_harmonics = n_harmonics,
    harmonic_amps = harmonic_amps, harmonic_offsets = harmonic_offsets,
    harmonic_periods = harmonic_periods,
    annual_jitter_sd = annual_jitter_sd, daily_noise_sd = daily_noise_sd,
    daily_noise_clamp = config$daily_noise_clamp,
    ss_event_prob = ss_event_prob, ss_boost = ss_boost,
    ss_shape = ss_shape, ss_scale = ss_scale,
    has_intervention = has_intervention,
    interv_mode = policy$mode,
    iso_on_day = start_day, iso_off_day = end_day, iso_active = iso_active,
    quar_on_day = start_day, quar_off_day = end_day, quar_active = quar_active,
    pep_on_day = start_day, pep_off_day = end_day, pep_active = pep_active,
    prep_on_day = start_day, prep_off_day = end_day,
    prep_start_rate = prep_start, prep_active = prep_eff,
    school_closure = school_red, work_closure = work_red,
    on_threshold = on_threshold, off_threshold = off_threshold,
    trigger_delay = trigger_delay, contact_reduction = contact_reduction,
    endemic = structure$has_demography,
    birth_rate = birth_rate, death_rate = death_rate,
    importation_rate = importation_rate,
    seed_frac = seed_frac,
    seed_mode = seed_mode, n_seed = n_seed, seed_comp = seed_comp, novel = population$novel,
    immune_frac = population$immune_frac,
    hosp_prob = hosp_prob, death_prob = death_prob,
    detection = detection, report_delay = report_delay, r_nb = r_nb
  )
}


################################################################################
# Mechanism-specific samplers
################################################################################

.finish_params <- function(base_list, epi, struct, scenario, structure) {
  age_hosp  <- pmin(0.95, scenario$hosp_prob  * struct$age_multipliers)
  age_death <- pmin(0.95, scenario$death_prob * struct$age_multipliers)
  c(base_list, epi,
    struct[c("P", "A", "C", "Mp", "pop_fracs", "age_multipliers")],
    scenario,
    list(age_hosp_probs = age_hosp, age_death_probs = age_death,
         structure = structure))
}


sample_params_base <- function(structure, population, policy, config) {
  p  <- config$r0_base
  R0 <- .draw_r0(p)

  epi <- sample_shared_epi(structure, config)
  struct <- sample_population_structure(population)

  if (structure$has_recovery) {
    ip <- compute_ip(epi$infectious_days, epi$latent_days, epi$p_asym,
                       epi$asymp_trans, epi$kappa, structure$has_latent)
    beta_raw <- R0 / max(ip * struct$structural_mult, 1e-8)
    beta <- pmax(config$beta_bounds[1], pmin(beta_raw, config$beta_bounds[2]))
    beta_clipped  <- !isTRUE(all.equal(beta, beta_raw))
    R0_realized   <- beta * ip * struct$structural_mult
  } else {
    #messages when beta is not defined
    if (is.null(config$beta_direct))
      stop("`has_recovery = FALSE` (SI model) requires `beta_direct` in ",
           "reservoir_config(): the R0 -> beta inversion diverges when gamma = 0.")
    ip <- Inf; beta <- config$beta_direct; R0 <- NA_real_
    beta_clipped <- FALSE; R0_realized <- NA_real_
  }

  scenario <- sample_scenario(beta, structure, population, policy, config)
  .finish_params(list(mechanism = "base", R0 = R0, ip = ip, beta = beta,
                      beta_clipped = beta_clipped, R0_realized = R0_realized,
                      structural_mult = struct$structural_mult),
                 epi, struct, scenario, structure)
}


sample_params_vector <- function(structure, population, policy, config) {
  p  <- config$r0_vector
  R0 <- .draw_r0(p)
  v  <- config$vector

  epi <- sample_shared_epi(structure, config)
  struct <- sample_population_structure(population)

  sigma_v <- .logunif(v$sigma_v_logrange)
  mu_v  <- .logunif(v$mu_v_logrange)
  a_bite <- .runif_r(v$a_bite_range)
  b_h <- .runif_r(v$b_h_range)
  b_v <- .runif_r(v$b_v_range)
  Nv_init_frac <- .runif_r(v$Nv_init_frac_range)
  Iv_frac <- .runif_r(v$Iv_frac_range)

  # Age-specific relative biting exposure
  A_n <- struct$A
  w_raw <- exp(rnorm(A_n, 0, v$bite_age_sdlog %||% 0.3))
  age_fracs <- struct$pop_fracs %||% rep(1 / A_n, A_n)
  bite_age_weights <- w_raw / sum(w_raw * age_fracs)

  if (structure$has_recovery) {
    ip <- compute_ip(epi$infectious_days, epi$latent_days, epi$p_asym,
                     epi$asymp_trans, epi$kappa, structure$has_latent)
  
    time_Iv <- (sigma_v / max(sigma_v + mu_v, 1e-8)) * (1 / max(mu_v, 1e-8))
    k2 <- a_bite * a_bite * b_h * b_v * time_Iv * ip * Nv_init_frac
    beta_raw <- R0 / max(sqrt(k2), 1e-8)
    beta <- pmax(config$beta_bounds[1], pmin(beta_raw, config$beta_bounds[2]))
    beta_clipped  <- !isTRUE(all.equal(beta, beta_raw))
    R0_realized   <- beta * sqrt(k2)

    # All transmission is vector-mediated.
    R0_h2h <- 0
    R0_vector <- R0
  } else {
    if (is.null(config$beta_direct))
      stop("`has_recovery = FALSE` requires `beta_direct` in reservoir_config()")
    ip <- Inf; beta <- config$beta_direct; R0 <- NA_real_
    beta_clipped <- FALSE; R0_realized <- NA_real_
    R0_h2h <- NA_real_; R0_vector <- NA_real_
  }

  scenario <- sample_scenario(beta, structure, population, policy, config)
  .finish_params(list(mechanism = "vector", R0 = R0, ip = ip, beta = beta,
                      beta_clipped = beta_clipped, R0_realized = R0_realized,
                      R0_h2h = R0_h2h, R0_vector = R0_vector,
                      structural_mult = struct$structural_mult,
                      sigma_v = sigma_v, mu_v = mu_v, a_bite = a_bite,
                      b_h = b_h, b_v = b_v,
                      Nv_init_frac = Nv_init_frac, Iv_frac = Iv_frac,
                      bite_age_weights = bite_age_weights),
                 epi, struct, scenario, structure)
}


sample_params_waterborne <- function(structure, population, policy, config) {
  p <- config$r0_waterborne
  R0 <- .draw_r0(p)
  w <- config$waterborne

  waterborne_frac <- .runif_r(w$waterborne_frac_range)
  contact_frac <- 1 - waterborne_frac

  epi <- sample_shared_epi(structure, config)
  struct <- sample_population_structure(population)

  eta_I <- .logunif(w$eta_I_logrange)
  eta_A <- .logunif(w$eta_A_logrange)
  alpha_env <- .runif_r(w$alpha_env_range)
  mu_w <- .logunif(w$mu_w_logrange)

  # delta_eff depends on N, so look at pop+soize first
  N_pop <- if (!is.null(population$pop_size)) as.numeric(population$pop_size)
           else round(exp(runif(1, log(population$pop_size_range[1]),
                                   log(population$pop_size_range[2]))))

  if (structure$has_recovery) {
    ip <- compute_ip(epi$infectious_days, epi$latent_days, epi$p_asym,
                     epi$asymp_trans, epi$kappa, structure$has_latent)
    R0_contact <- R0 * contact_frac
    beta_raw <- R0_contact / max(ip * struct$structural_mult, 1e-8)
    beta_contact <- pmax(config$beta_bounds[1],
                         pmin(beta_raw, config$beta_bounds[2]))
    # Either bound bi
    beta_clipped <- !isTRUE(all.equal(beta_contact, beta_raw))

    R0_water <- R0 * waterborne_frac
    shed_rate <- eta_I * (1 - epi$p_asym) + eta_A * alpha_env * epi$p_asym
    # W accumulates ABSOLUTE shedding,
    mean_W <- shed_rate * ip / max(mu_w, 1e-8)
    delta_raw <- R0_water / max(mean_W * N_pop, 1e-8)
    delta_eff <- pmax(w$delta_eff_bounds[1],
                      pmin(delta_raw, w$delta_eff_bounds[2]))
    delta_clipped <- !isTRUE(all.equal(delta_eff, delta_raw))
    R0_realized   <- beta_contact * ip * struct$structural_mult +
                     delta_eff * mean_W * N_pop
  } else {
    if (is.null(config$beta_direct))
      stop("`has_recovery = FALSE` requires `beta_direct` in reservoir_config()")
    ip <- Inf; beta_contact <- config$beta_direct
    R0 <- NA_real_; R0_contact <- NA_real_; R0_water <- NA_real_
    delta_eff <- w$delta_eff_bounds[1]
    beta_clipped <- FALSE; delta_clipped <- FALSE; R0_realized <- NA_real_
  }

  population$pop_size <- N_pop   # pin, so scenario uses the same N as delta_eff
  scenario <- sample_scenario(beta_contact, structure, population, policy, config)
  .finish_params(list(mechanism = "waterborne", R0 = R0,
                      beta_clipped = beta_clipped, delta_clipped = delta_clipped,
                      R0_realized = R0_realized,
                      R0_contact = R0_contact, R0_water = R0_water,
                      waterborne_frac = waterborne_frac,
                      ip = ip, beta = beta_contact,
                      structural_mult = struct$structural_mult,
                      eta_I = eta_I, eta_A = eta_A, w_unit = w$w_unit %||% 1.0,
                      alpha_env = alpha_env, mu_w = mu_w, delta_eff = delta_eff),
                 epi, struct, scenario, structure)
}


################################################################################
# sample_params — unified entry point
################################################################################
#' Draw a hierarchical parameter set for one run
#'
#' @param mechanism "base", "vector", or "waterborne".
#' @param structure A `model_structure()` object, or an already-realised
#'   structure list from `.sample_structure()`.
#' @param population A `population_spec()`.
#' @param policy An `intervention_policy()`.
#' @param config A `reservoir_config()`.
#'
#' @export
sample_params <- function(mechanism  = "base",
                          structure  = model_structure(),
                          population = population_spec(),
                          policy     = intervention_policy(),
                          config     = reservoir_config()) {

  if (inherits(structure, "reservoir_structure"))
    structure <- .sample_structure(structure)

  switch(mechanism,
    base  = sample_params_base(structure, population, policy, config),
    vector  = sample_params_vector(structure, population, policy, config),
    waterborne = sample_params_waterborne(structure, population, policy, config),
    stop("mechanism must be 'base', 'vector', or 'waterborne'")
  )
}


################################################################################
# DIAGNOSTICS
################################################################################
#' give overview of distributions from sampler
#'
#' @export
check_sampler <- function(mechanism = "base", n = 500,
                          structure  = model_structure(),
                          population = population_spec(),
                          policy     = intervention_policy(),
                          config     = reservoir_config()) {

  cat(sprintf("\n=== %s sampler (n=%d) ===\n", toupper(mechanism), n))
  samples <- lapply(seq_len(n), function(i)
    sample_params(mechanism, structure, population, policy, config))

  pct <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return("(all non-finite)")
    sprintf("%.2f  [%.2f, %.2f]", mean(x), quantile(x, 0.05), quantile(x, 0.95))
  }
  g <- function(f) sapply(samples, f)

  cat(sprintf("  R0:              %s\n", pct(g(function(s) s$R0))))
  cat(sprintf("  beta:            %s\n", pct(g(function(s) s$beta))))
  cat(sprintf("  infectious_days: %s\n", pct(g(function(s) s$infectious_days))))
  cat(sprintf("  latent_days:     %s\n", pct(g(function(s) s$latent_days))))
  cat(sprintf("  p_asym:          %s\n", pct(g(function(s) s$p_asym))))
  Ns <- g(function(s) s$N)
  cat(sprintf("  N (pop):         mean %.0f  [%.0f, %.0f]\n",
              mean(Ns), quantile(Ns, 0.05), quantile(Ns, 0.95)))

  sk <- table(g(function(s) s$structure$skeleton))
  cat("  skeletons:\n")
  for (nm in names(sk))
    cat(sprintf("    %-10s %d/%d (%.0f%%)\n", nm, sk[[nm]], n, 100 * sk[[nm]] / n))

  invisible(samples)
}
