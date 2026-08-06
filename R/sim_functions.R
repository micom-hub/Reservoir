################################################################################
# sim_functions.R
#
#Rt calculations
################################################################################


################################################################################
# run_sim - simulator + post-processing + Rt in one call
################################################################################
#' Run one simulation with post-processing and Rt
#'
#' @export
run_sim <- function(cfg,
                    hosp_prob  = 0.05,
                    death_prob = 0.01,
                    seed       = NULL,
                    mechanism  = "base",
                    config     = reservoir_config()) {

  if (!is.null(seed)) set.seed(seed)

  result <- .dispatch_engine(mechanism, cfg)

  out <- run_postprocessing(
    result   = result, params = cfg$params,
    P = cfg$P, A = cfg$A, num_days = cfg$tf_days,
    age_hosp_probs = hosp_prob, age_death_probs = death_prob,
    config = config
  )

  Rt <- compute_daily_rt(result, cfg, mechanism = mechanism)

  total_cases <- sum(out$daily_cases_true)
  nz_Rt <- Rt[is.finite(Rt) & Rt > 0]

  c(out, list(
    daily_Rt    = Rt,
    final_state = result$final_state,
    cumul       = result$cumul_infections,
    summary = list(
      total_cases  = total_cases,
      attack_rate  = round(total_cases / max(cfg$N, 1), 3),
      peak_day     = which.max(out$daily_cases_true),
      peak_cases   = max(out$daily_cases_true),
      total_hosp   = sum(out$daily_hosp_true),
      total_deaths = sum(out$daily_death_true),
      peak_Rt      = if (length(nz_Rt)) round(max(nz_Rt), 2) else NA,
      final_Rt     = if (length(nz_Rt)) round(tail(nz_Rt, 1), 2) else NA
    )
  ))
}


################################################################################
# .dispatch_engine - pick the compiled engine for a mechanism
################################################################################
.dispatch_engine <- function(mechanism, cfg) {
  switch(mechanism,
    base = simulate_structured_atl_cpp(cfg),
    vector = simulate_vectorborne_atl_cpp(cfg),
    waterborne = simulate_env_atl_cpp(cfg),
    stop("unknown mechanism: ", mechanism)
  )
}


################################################################################
# .as_sim - normalize the two run shapes into one
#
# run_sim() returns a FLAT object (daily_cases_true, daily_hosp_true, ...).
# generate_dataset() returns a run with `daily_true` and the observation
# streams nested under `noise_full`. Users naturally pass either to plot_sim()
# and summary_sim(), so both are accepted here.
################################################################################
.as_sim <- function(x) {
  if (!is.null(x$daily_cases_true)) return(x)          # already flat
  if (is.null(x$daily_true))
    stop("object does not look like a Reservoir run: no `daily_cases_true` ",
         "or `daily_true` field")

  out <- list(daily_cases_true = x$daily_true,
              daily_Rt         = x$daily_Rt,
              summary          = x$summary)

  # cheap noise stream is always present
  if (!is.null(x$noise_limited$daily_cases))
    out$daily_cases_reported <- x$noise_limited$daily_cases

  # full observation bundle, when the run scope produced one
  nf <- x$noise_full
  if (!is.null(nf)) {
    for (nm in c("daily_cases_reported", "daily_hosp_true", "daily_hosp_reported",
                 "daily_death_true", "daily_death_reported", "hosp_occupancy",
                 "wastewater_raw", "wastewater_normalized", "syndromic_signal"))
      if (!is.null(nf[[nm]])) out[[nm]] <- nf[[nm]]
  }
  out
}


################################################################################
# summary_sim / plot_sim
################################################################################
#' @export
summary_sim <- function(sim_out) {
  sim_out <- .as_sim(sim_out)
  s <- sim_out$summary
  cat(sprintf("Attack rate:   %.1f%%\n",            s$attack_rate * 100))
  cat(sprintf("Peak day:      day %d (%d cases)\n", s$peak_day, s$peak_cases))
  cat(sprintf("Total hosp:    %d\n",                s$total_hosp))
  cat(sprintf("Total deaths:  %d\n",                s$total_deaths))
  cat(sprintf("Peak Rt:       %.2f\n",              s$peak_Rt))
  cat(sprintf("Final Rt:      %.2f\n",              s$final_Rt))
  invisible(s)
}


#' @export
plot_sim <- function(sim_out, what = "all", title = NULL) {
  sim_out <- .as_sim(sim_out)
  valid <- c("all", "cases", "hosp", "deaths", "rt", "wastewater", "syndromic")
  if (!what %in% valid) stop("what must be one of: ", paste(valid, collapse = ", "))

  plots <- list(
    cases = function() {
      plot(sim_out$daily_cases_true, type = "l", lwd = 2,
           xlab = "Day", ylab = "Daily cases", main = "Epidemic curve")
      if (!is.null(sim_out$daily_cases_reported))
        lines(sim_out$daily_cases_reported, col = "steelblue", lty = 2, lwd = 1.5)
      legend("topright", c("True", "Reported"), col = c("black", "steelblue"),
             lty = 1:2, lwd = 1.5, bty = "n", cex = 0.8)
    },
    hosp = function() {
      ymax <- max(c(sim_out$daily_hosp_true, sim_out$hosp_occupancy), na.rm = TRUE)
      plot(sim_out$daily_hosp_true, type = "l", lwd = 2, col = "darkred",
           ylim = c(0, ymax), xlab = "Day", ylab = "Count", main = "Hospitalisations")
      lines(sim_out$daily_hosp_reported, col = "salmon",  lty = 2, lwd = 1.5)
      lines(sim_out$hosp_occupancy,      col = "darkred", lty = 3, lwd = 1.5)
      legend("topright", c("Admissions (true)", "Admissions (rep.)", "Occupancy"),
             col = c("darkred", "salmon", "darkred"), lty = 1:3, lwd = 1.5,
             bty = "n", cex = 0.8)
    },
    deaths = function() {
      plot(sim_out$daily_death_true, type = "l", lwd = 2,
           xlab = "Day", ylab = "Deaths", main = "Deaths")
      lines(sim_out$daily_death_reported, col = "gray50", lty = 2, lwd = 1.5)
    },
    rt = function() {
      ymax <- max(sim_out$daily_Rt, 2, na.rm = TRUE)
      plot(sim_out$daily_Rt, type = "l", lwd = 2, col = "purple",
           ylim = c(0, ymax), xlab = "Day", ylab = "Rt",
           main = "Effective reproduction number")
      abline(h = 1, lty = 2, col = "red")
    },
    wastewater = function() {
      plot(sim_out$wastewater_normalized, type = "l", lwd = 2, col = "darkgreen",
           xlab = "Day", ylab = "Normalised concentration", main = "Wastewater signal")
    },
    syndromic = function() {
      plot(sim_out$syndromic_signal, type = "l", lwd = 2, col = "orange",
           xlab = "Day", ylab = "Consultations", main = "Syndromic surveillance")
    }
  )

  # drop panels whose stream this run does not carry
  avail <- c(cases = !is.null(sim_out$daily_cases_true),
             hosp  = !is.null(sim_out$daily_hosp_true),
             deaths= !is.null(sim_out$daily_death_true),
             rt    = !is.null(sim_out$daily_Rt) && any(is.finite(sim_out$daily_Rt)),
             wastewater = !is.null(sim_out$wastewater_normalized),
             syndromic  = !is.null(sim_out$syndromic_signal))
  if (what != "all" && !isTRUE(avail[[what]]))
    stop("this run has no '", what, "' stream (was it generated with a ",
         "narrower `scope`?)")

  if (what == "all") {
    plots <- plots[names(plots) %in% names(avail)[avail]]
    if (!length(plots)) stop("nothing to plot")
    nr <- if (length(plots) <= 3) 1 else 2
    nc <- ceiling(length(plots) / nr)
    old <- par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1),
               oma = if (!is.null(title)) c(0, 0, 2, 0) else c(0, 0, 0, 0))
    on.exit(par(old))
    for (fn in plots) try(fn(), silent = TRUE)
    if (!is.null(title)) mtext(title, outer = TRUE, cex = 1.1, font = 2)
  } else {
    plots[[what]]()
    if (!is.null(title)) title(main = title)
  }
  invisible(sim_out)
}


#' @export
run_replicates <- function(cfg, n_reps = 10L, hosp_prob = 0.05,
                           death_prob = 0.01, seed = NULL, mechanism = "base") {
  lapply(seq_len(n_reps), function(i) {
    s <- if (!is.null(seed)) seed + i else NULL
    run_sim(cfg, hosp_prob, death_prob, seed = s, mechanism = mechanism)
  })
}


################################################################################
# Structured Rt helpers
################################################################################

.rt_driver <- function(cfg, name, day, default = 0.0) {
  v <- cfg[[name]]
  if (is.null(v)) return(default)
  v[min(day, length(v))]
}

.rt_block_slice <- function(state, block, offset) {
  state[(offset * block + 1L):((offset + 1L) * block)]
}

.rt_extract_human_state <- function(state, block) {
  out <- list(
    S     = .rt_block_slice(state, block, 0L),
    Sprep = .rt_block_slice(state, block, 1L),
    Spost = .rt_block_slice(state, block, 2L),
    Qs    = .rt_block_slice(state, block, 3L),
    L     = .rt_block_slice(state, block, 4L),
    QL    = .rt_block_slice(state, block, 5L),
    Is    = .rt_block_slice(state, block, 6L),
    Ia    = .rt_block_slice(state, block, 7L),
    Iso   = .rt_block_slice(state, block, 8L),
    R     = .rt_block_slice(state, block, 9L),
    R2    = .rt_block_slice(state, block, 10L),
    R3    = .rt_block_slice(state, block, 11L)
  )
  out$N <- Reduce(`+`, out)
  out
}

.rt_cell_idx <- function(p, a, A) (p - 1L) * A + a

# Index of infected-subsystem state `comp` (1-based within the subsystem)
# for cell `cell`. Layout is [state1 cells..., state2 cells..., ...].
.rt_idx <- function(comp, cell, block) (comp - 1L) * block + cell


# -----------------------------------------------------------------------------
# .rt_subsystem — which infected states exist, and where infections land
#
# Returns:
#   n_states   number of infected states per cell
#   entry      integer vector of subsystem indices that new infections enter
#   entry_w    weights for those entries (must sum to 1)
#   is_idx     subsystem index of Is (transmits at weight 1)
#   ia_idx     subsystem index of Ia (transmits at asymp_trans), or NA
# -----------------------------------------------------------------------------
.rt_subsystem <- function(has_latent, p_asym) {
  if (has_latent) {
    # {L, Is, Ia, Iso, QL}
    list(n_states = 5L, entry = 1L, entry_w = 1.0, is_idx = 2L, ia_idx = 3L)
  } else {
    # {Is, Ia} — infections branch at entry
    list(n_states = 2L, entry = c(1L, 2L),
         entry_w = c(1 - p_asym, p_asym), is_idx = 1L, ia_idx = 2L)
  }
}


################################################################################
# .rt_V — transition matrix of the infected subsystem
################################################################################
.rt_V <- function(block, sub, sigma, gamma, kappa, p_asym,
                  iso_t, quar_t, pep_t, death_rate, extra_dim = 0L) {

  dim <- sub$n_states * block + extra_dim
  V <- matrix(0, dim, dim)

  if (sub$n_states == 5L) {
    mu_L   <- sigma + quar_t + pep_t + death_rate
    mu_Is  <- gamma + iso_t  + death_rate
    mu_Ia  <- kappa + iso_t  + gamma + death_rate
    mu_Iso <- gamma + death_rate
    mu_QL  <- sigma + pep_t  + death_rate

    for (cell in seq_len(block)) {
      iL   <- .rt_idx(1L, cell, block); iIs  <- .rt_idx(2L, cell, block)
      iIa  <- .rt_idx(3L, cell, block); iIso <- .rt_idx(4L, cell, block)
      iQL  <- .rt_idx(5L, cell, block)

      V[iL, iL]     <-  mu_L
      V[iIs, iIs]   <-  mu_Is
      V[iIa, iIa]   <-  mu_Ia
      V[iIso, iIso] <-  mu_Iso
      V[iQL, iQL]   <-  mu_QL

      V[iIs,  iL]  <- -sigma * (1 - p_asym)
      V[iIa,  iL]  <- -sigma * p_asym
      V[iQL,  iL]  <- -quar_t
      V[iIs,  iIa] <- -kappa
      V[iIso, iIs] <- -iso_t
      V[iIso, iIa] <- -iso_t
      V[iIa,  iQL] <- -sigma * p_asym
      V[iIso, iQL] <- -sigma * (1 - p_asym)
    }
  } else {
    # {Is, Ia}: no latent stage, no QL, and Iso is a pure sink (does not
    # transmit and does not feed back), so it is excluded from the subsystem.
    mu_Is <- gamma + iso_t + death_rate
    mu_Ia <- kappa + iso_t + gamma + death_rate

    for (cell in seq_len(block)) {
      iIs <- .rt_idx(1L, cell, block); iIa <- .rt_idx(2L, cell, block)
      V[iIs, iIs] <-  mu_Is
      V[iIa, iIa] <-  mu_Ia
      V[iIs, iIa] <- -kappa          # Ia -> Is
    }
  }
  V
}


################################################################################
# .rt_add_contact_F — human-to-human new-infection matrix
#
# Deposits new infections into sub$entry with weights sub$entry_w, which is
# what makes this work for both the latent and no-latent subsystems.
################################################################################
.rt_add_contact_F <- function(F_mat, h, cfg, sub, beta_t, s_eff, scale = 1.0) {
  P <- cfg$P; A <- cfg$A
  block <- P * A
  C  <- cfg$C; Mp <- cfg$Mp
  asymp_trans <- cfg$params$asymp_trans

  for (p in seq_len(P)) for (a in seq_len(A)) {
    target_cell <- .rt_cell_idx(p, a, A)

    for (q in seq_len(P)) {
      mpq <- Mp[p, q]
      if (mpq == 0) next
      for (b in seq_len(A)) {
        source_cell <- .rt_cell_idx(q, b, A)
        N_source <- h$N[source_cell]
        if (N_source <= 0) next

        contact <- beta_t * scale * mpq * C[a, b] * s_eff[target_cell] / N_source

        src_is <- .rt_idx(sub$is_idx, source_cell, block)
        src_ia <- if (!is.na(sub$ia_idx)) .rt_idx(sub$ia_idx, source_cell, block) else NA

        for (k in seq_along(sub$entry)) {
          tgt <- .rt_idx(sub$entry[k], target_cell, block)
          w   <- sub$entry_w[k]
          F_mat[tgt, src_is] <- F_mat[tgt, src_is] + w * contact
          if (!is.na(src_ia))
            F_mat[tgt, src_ia] <- F_mat[tgt, src_ia] + w * asymp_trans * contact
        }
      }
    }
  }
  F_mat
}


.rt_spectral_radius <- function(F_mat, V) {
  V_inv <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(V_inv)) return(NA_real_)
  eig <- tryCatch(eigen(F_mat %*% V_inv, only.values = TRUE)$values,
                  error = function(e) NULL)
  if (is.null(eig)) return(NA_real_)
  max(Re(eig))
}


.rt_common_day <- function(result, cfg, day) {
  block <- cfg$P * cfg$A
  state <- result$state_over_time[day, ]
  h <- .rt_extract_human_state(state, block)
  prep_t <- .rt_driver(cfg, "prep_eff_daily", day, default = 0.0)

  # Reactive runs gate interventions inside the engine; mirror that here so
  # Rt reflects what actually happened rather than the supplied driver value.
  active <- 1.0
  if (isTRUE(cfg$reactive_mode == 1L) && !is.null(result$intervention_active))
    active <- result$intervention_active[min(day, length(result$intervention_active))]

  beta_t <- .rt_driver(cfg, "beta_eff_daily", day)
  if (isTRUE(cfg$reactive_mode == 1L) && active > 0)
    beta_t <- beta_t * (cfg$reactive_contact_mult %||% 1.0)

  s_eff <- h$S + (1 - prep_t * active) * h$Sprep

  list(
    block = block, state = state, human = h, s_eff = s_eff,
    beta_t = beta_t,
    iso_t  = .rt_driver(cfg, "iso_rate_daily",  day) * active,
    quar_t = .rt_driver(cfg, "quar_rate_daily", day) * active,
    pep_t  = .rt_driver(cfg, "pep_rate_daily",  day) * active
  )
}


################################################################################
# compute_daily_rt — dispatches on mechanism and structure
#
# Correct for SIS with no change: omega and R never enter F or V.
################################################################################
#' Daily effective reproduction number via the next-generation matrix
#'
#' @export
compute_daily_rt <- function(result, cfg, mechanism = "base") {
  has_latent <- !isTRUE(cfg$no_latent == 1L)
  sub <- .rt_subsystem(has_latent, cfg$params$p_asym)

  switch(mechanism,
    base       = .rt_loop_base(result, cfg, sub),
    vector     = .rt_loop_vector(result, cfg, sub),
    waterborne = .rt_loop_waterborne(result, cfg, sub),
    stop("unknown mechanism: ", mechanism)
  )
}

# Back-compatible named entry points
#' @export
compute_daily_rt_vector <- function(result, cfg)
  compute_daily_rt(result, cfg, mechanism = "vector")

#' @export
compute_daily_rt_waterborne <- function(result, cfg)
  compute_daily_rt(result, cfg, mechanism = "waterborne")


.rt_loop_base <- function(result, cfg, sub) {
  Rt <- numeric(cfg$tf_days)
  for (day in seq_len(cfg$tf_days)) {
    d <- .rt_common_day(result, cfg, day)
    dim_h <- sub$n_states * d$block
    F_mat <- matrix(0, dim_h, dim_h)
    F_mat <- .rt_add_contact_F(F_mat, d$human, cfg, sub, d$beta_t, d$s_eff)
    V <- .rt_V(d$block, sub, cfg$params$sigma, cfg$params$gamma,
               cfg$params$kappa, cfg$params$p_asym,
               d$iso_t, d$quar_t, d$pep_t, cfg$params$death_rate)
    Rt[day] <- .rt_spectral_radius(F_mat, V)
  }
  Rt
}


.rt_loop_vector <- function(result, cfg, sub) {
  block   <- cfg$P * cfg$A
  M_human <- 13L * block
  Rt <- numeric(cfg$tf_days)

  for (day in seq_len(cfg$tf_days)) {
    d <- .rt_common_day(result, cfg, day)
    dim_h <- sub$n_states * d$block
    iEv <- dim_h + 1L; iIv <- dim_h + 2L
    dim_total <- dim_h + 2L

    seasonal_t <- .rt_driver(cfg, "seasonal_daily", day, default = 1.0)
    f_interv_t <- .rt_driver(cfg, "f_interv_daily", day, default = 1.0)
    beta_route <- d$beta_t * seasonal_t * f_interv_t

    F_mat <- matrix(0, dim_total, dim_total)
    F_mat <- .rt_add_contact_F(F_mat, d$human, cfg, sub, d$beta_t, d$s_eff,
                               scale = seasonal_t * f_interv_t)

    Sv <- d$state[M_human + 1L]; Ev <- d$state[M_human + 2L]; Iv <- d$state[M_human + 3L]
    Nv <- Sv + Ev + Iv
    Nh_total <- sum(d$human$N)

    for (cell in seq_len(d$block)) {
      # vector -> human: enters wherever new infections enter
      for (k in seq_along(sub$entry)) {
        tgt <- .rt_idx(sub$entry[k], cell, d$block)
        F_mat[tgt, iIv] <- F_mat[tgt, iIv] +
          sub$entry_w[k] * beta_route * cfg$params$a_bite * cfg$params$b_h *
          d$s_eff[cell] / max(Nv, 1)
      }
      # human -> vector
      src_is <- .rt_idx(sub$is_idx, cell, d$block)
      F_mat[iEv, src_is] <- beta_route * cfg$params$a_bite * cfg$params$b_v *
        Sv / max(Nh_total, 1)
      if (!is.na(sub$ia_idx)) {
        src_ia <- .rt_idx(sub$ia_idx, cell, d$block)
        F_mat[iEv, src_ia] <- beta_route * cfg$params$asymp_trans *
          cfg$params$a_bite * cfg$params$b_v * Sv / max(Nh_total, 1)
      }
    }

    V <- .rt_V(d$block, sub, cfg$params$sigma, cfg$params$gamma,
               cfg$params$kappa, cfg$params$p_asym,
               d$iso_t, d$quar_t, d$pep_t, cfg$params$death_rate,
               extra_dim = 2L)
    V[iEv, iEv] <- cfg$params$sigma_v + cfg$params$mu_v
    V[iIv, iIv] <- cfg$params$mu_v
    V[iIv, iEv] <- -cfg$params$sigma_v

    Rt[day] <- .rt_spectral_radius(F_mat, V)
  }
  Rt
}


.rt_loop_waterborne <- function(result, cfg, sub) {
  Rt <- numeric(cfg$tf_days)
  for (day in seq_len(cfg$tf_days)) {
    d <- .rt_common_day(result, cfg, day)
    dim_h <- sub$n_states * d$block

    f_contact_t <- .rt_driver(cfg, "f_contact_daily", day, default = 1.0)
    f_water_t   <- .rt_driver(cfg, "f_water_daily",   day, default = 1.0)
    delta_t     <- .rt_driver(cfg, "delta_eff_daily", day)

    F_mat <- matrix(0, dim_h, dim_h)
    F_mat <- .rt_add_contact_F(F_mat, d$human, cfg, sub, d$beta_t, d$s_eff,
                               scale = f_contact_t)

    # W collapsed at quasi-steady state
    w_per_Is <- cfg$params$eta_I / max(cfg$params$mu_w, 1e-8)
    w_per_Ia <- cfg$params$eta_A * cfg$params$alpha_env / max(cfg$params$mu_w, 1e-8)

    for (target_cell in seq_len(d$block)) {
      water_Is <- delta_t * f_water_t * d$s_eff[target_cell] * w_per_Is
      water_Ia <- delta_t * f_water_t * d$s_eff[target_cell] * w_per_Ia
      for (source_cell in seq_len(d$block)) {
        src_is <- .rt_idx(sub$is_idx, source_cell, d$block)
        src_ia <- if (!is.na(sub$ia_idx)) .rt_idx(sub$ia_idx, source_cell, d$block) else NA
        for (k in seq_along(sub$entry)) {
          tgt <- .rt_idx(sub$entry[k], target_cell, d$block)
          F_mat[tgt, src_is] <- F_mat[tgt, src_is] + sub$entry_w[k] * water_Is
          if (!is.na(src_ia))
            F_mat[tgt, src_ia] <- F_mat[tgt, src_ia] + sub$entry_w[k] * water_Ia
        }
      }
    }

    V <- .rt_V(d$block, sub, cfg$params$sigma, cfg$params$gamma,
               cfg$params$kappa, cfg$params$p_asym,
               d$iso_t, d$quar_t, d$pep_t, cfg$params$death_rate)
    Rt[day] <- .rt_spectral_radius(F_mat, V)
  }
  Rt
}
