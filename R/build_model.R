################################################################################
# build_model.R
#
# This funciton does the following:
#   1. Takes a compartment + reaction specification
#   2. Calls generate_slots.R to build each C++ code fragment
#   3. Fills those fragments into .cpp files
#   4. Writes the result to a temp file
#   5. Calls sourceCpp() to compile and load it
#   6. Returns a custom setup_sim() function pre-configured for your model
#
################################################################################



################################################################################
# validate_spec
#
# Checks the compartment/reaction specification for common mistakes
# before attempting code generation.
################################################################################
validate_spec <- function(compartments, infectious, in_denominator, reactions) {

  # cumul must be present, needed for daily case counting
  if (!"cumul" %in% compartments)
    stop("compartments must include 'cumul' (cumulative infection tracker)")

  # All infectious compartments must exist
  inf_comps <- vapply(infectious, `[[`, character(1), "comp")
  bad_inf <- setdiff(inf_comps, compartments)
  if (length(bad_inf) > 0)
    stop("infectious compartments not in compartments list: ",
         paste(bad_inf, collapse=", "))

  # All denominator compartments must exist
  bad_den <- setdiff(in_denominator, compartments)
  if (length(bad_den) > 0)
    stop("in_denominator compartments not in compartments list: ",
         paste(bad_den, collapse=", "))

  # cumul should not be in the denominator
  if ("cumul" %in% in_denominator)
    stop("'cumul' should not be in in_denominator — it is a counter, not a population")

  # All reaction from/to compartments must exist
  for (i in seq_along(reactions)) {
    r <- reactions[[i]]
    type <- if (!is.null(r$type)) r$type else "flow"
    if (type %in% c("birth", "death")) next

    if (is.null(r$from))
      stop("Reaction ", i, " has no 'from' compartment")

    if (!r$from %in% compartments)
      stop("Reaction ", i, " 'from' compartment '", r$from,
           "' not in compartments list")

    if (type == "flow") {
      if (is.null(r$to))
        stop("Reaction ", i, " has no 'to' compartment")
      if (!r$to %in% compartments)
        stop("Reaction ", i, " 'to' compartment '", r$to,
             "' not in compartments list")
      if (is.null(r$rate))
        stop("Reaction ", i, " (", r$from, " -> ", r$to, ") has no 'rate'")
    }

    if (type == "infection") {
      if (is.null(r$to))
        stop("Infection reaction ", i, " has no 'to' targets")
      bad_to <- setdiff(r$to, compartments)
      if (length(bad_to) > 0)
        stop("Infection reaction ", i, " 'to' targets not in compartments: ",
             paste(bad_to, collapse=", "))
    }
  }

  # Must have exactly one birth and one death reaction
  types <- vapply(reactions, function(r)
    if (!is.null(r$type)) r$type else "flow", character(1))
  if (sum(types == "birth") != 1)
    stop("Exactly one 'birth' reaction required (got ",
         sum(types == "birth"), ")")
  if (sum(types == "death") != 1)
    stop("Exactly one 'death' reaction required (got ",
         sum(types == "death"), ")")

  # Must have at least one infection reaction
  if (sum(types == "infection") < 1)
    stop("At least one 'infection' reaction required")

  invisible(TRUE)
}


################################################################################
# generate_cpp
#
# Generates all slot strings and fills the template
# Returns the complete C++ source as a character string
################################################################################
generate_cpp <- function(compartments, infectious, in_denominator,
                          reactions, template_path) {

  slots <- list(
    N_BLOCKS = gen_n_blocks(compartments),
    OFFSET_PARAMS = gen_offset_params(compartments),
    LOCAL_INDEX_VARS = gen_local_index_vars(compartments),
    REACTION_LINES = gen_reaction_lines(reactions, compartments),
    LAMBDA_VAR_DECLS = gen_lambda_var_decls(in_denominator, infectious),
    N_SUM = gen_n_sum(in_denominator),
    IEFF = gen_ieff(infectious),
    PROP_VAR_DECLS = gen_prop_var_decls(compartments),
    PROPENSITY_LINES = gen_propensity_lines(reactions, compartments),
    OFFSET_ASSIGNMENTS = gen_offset_assignments(compartments),
    OFFSET_ARGS = gen_offset_args(compartments),
    CUMUL_OFFSET_IDX = gen_cumul_offset_idx(compartments),
    DRIVER_READS = gen_driver_reads(reactions),
    DRIVER_DRV_LINES = gen_driver_drv_lines(reactions),
    PARAM_READS = gen_param_reads(reactions),
    RATE_PARAMS = gen_rate_params(reactions),
    RATE_ARGS = gen_rate_args(reactions)
  )

  fill_template(template_path, slots)
}


################################################################################
# make_setup_sim
#
# Returns a setup_sim() function pre-configured for the generated model.
# The returned function knows:
#   - How many blocks are in x0
#   - Which compartment offsets exist
#   - Which drivers are needed
#   - Which params are needed
################################################################################
make_setup_sim <- function(compartments, reactions) {

  n_blocks <- length(compartments)
  needed_drivers <- detect_needed_drivers(reactions)
  extra_params <- detect_extra_params(reactions)

  # Build offset lookup for x0 construction
  offsets <- setNames(seq_along(compartments) - 1L, compartments)

  function(
    N = 10000,
    P  = 1,
    A  = 1,
    tf_days = 365,
    beta = 0.3,
    infectious_days = 7,
    latent_days = 5,
    immunity_days = 365,
    p_asym = 0.3,
    asymp_trans = 0.5,
    kappa = 0.1,
    n_seed = 10,
    seed_comp = NULL,   # defaults to first non-cumul compartment
    birth_rate = 0.0,
    death_rate = 0.0,
    wave_change_days = NULL,
    wave_change_betas = NULL,
    use_seasonality = FALSE,
    seasonal_amplitude = 0.3,
    seasonal_peak_day  = 0,
    seasonal_period = 365,
    # Extra biological params (e.g. hosp_rate if you added hospitalisation)
    extra_param_values = list(),
    # Daily driver values scalar or length-tf_days vector
    iso_rate   = 0.0, iso_on_day   = NULL, iso_off_day   = NULL, iso_active   = 0.3,
    quar_rate  = 0.0, quar_on_day  = NULL, quar_off_day  = NULL, quar_active  = 0.2,
    pep_rate   = 0.0, pep_on_day   = NULL, pep_off_day   = NULL, pep_active   = 0.1,
    prep_eff   = 0.0, prep_on_day  = NULL, prep_off_day  = NULL, prep_active  = 0.7,
    leave_quar_rate = 0.1,
    leave_quar_prep_rate = 0.1,
    prep_start_rate = 0.0,
    prep_stop_rate = 0.0,
    C  = NULL,
    Mp = NULL,
    epsilon = 0.03,
    Ncritical = 10L,
    exactThreshold = 10.0,
    maxtau = Inf
  ) {

    block <- P * A

    # x
    if (is.null(seed_comp)) {
      seed_comp <- compartments[compartments != "cumul"][1]
    }
    if (!seed_comp %in% compartments)
      stop("seed_comp '", seed_comp, "' not in compartments")

    x0 <- numeric(n_blocks * block)
    x0[1] <- N - n_seed   # S (always first compartment)
    x0[offsets[seed_comp] * block + 1L] <- n_seed

    #beta_eff_daily
    if (!is.null(wave_change_days)) {
      if (length(wave_change_betas) != length(wave_change_days))
        stop("wave_change_days and wave_change_betas must be the same length")
      all_days  <- c(0L, as.integer(wave_change_days))
      all_betas <- c(beta, as.numeric(wave_change_betas))
      beta_vec  <- numeric(tf_days)
      for (day in seq_len(tf_days)) {
        seg           <- sum(all_days < day)
        beta_vec[day] <- all_betas[seg]
      }
    } else {
      beta_vec <- rep(beta, tf_days)
    }

    if (use_seasonality) {
      days         <- seq_len(tf_days) - 1L
      seasonal_vec <- pmax(0, 1 + seasonal_amplitude *
                             cos(2 * pi * (days - seasonal_peak_day) /
                                   seasonal_period))
    } else {
      seasonal_vec <- rep(1.0, tf_days)
    }
    beta_vec <- beta_vec * seasonal_vec

    # Intervention driver
    make_driver <- function(baseline, on_day, off_day, active) {
      vec <- rep(as.numeric(baseline), tf_days)
      if (!is.null(on_day)) {
        on_day  <- as.integer(on_day)
        off_day <- if (!is.null(off_day)) as.integer(off_day) else tf_days
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

    # Matrices for ages and populations
    if (is.null(C))  C  <- matrix(1.0, A, A)
    if (is.null(Mp)) Mp <- matrix(1.0 / P, P, P)

    # Parameter list building
    params <- list(
      sigma  = 1 / max(latent_days,     1e-8),
      gamma = 1 / max(infectious_days, 1e-8),
      omega  = 1 / max(immunity_days,   1e-8),
      p_asym = p_asym,
      asymp_trans = asymp_trans,
      kappa = kappa,
      birth_rate  = birth_rate,
      death_rate  = death_rate
    )
    # Append any extra params (e.g. hosp_rate, discharge_rate)
    for (nm in names(extra_param_values))
      params[[nm]] <- extra_param_values[[nm]]

    # Build cfg
    cfg <- list(
      P = P,
      A = A,
      tf_days = tf_days,
      C = C,
      Mp = Mp,
      x0 = x0,
      beta_eff_daily  = beta_vec,
      params = params,
      epsilon = epsilon,
      Ncritical = Ncritical,
      exactThreshold  = exactThreshold,
      maxtau = maxtau,
      N = N
    )

    # Add only the drivers actually needed by this model
    driver_map <- list(
      iso_rate_day  = make_driver(iso_rate,  iso_on_day,  iso_off_day,  iso_active),
      quar_rate_day = make_driver(quar_rate, quar_on_day, quar_off_day, quar_active),
      leave_quar_day = expand(leave_quar_rate,      "leave_quar_rate"),
      leave_quar_prep_day = expand(leave_quar_prep_rate, "leave_quar_prep_rate"),
      prep_start_rate_day = expand(prep_start_rate,      "prep_start_rate"),
      prep_stop_rate_day  = expand(prep_stop_rate,       "prep_stop_rate"),
      pep_rate_day = make_driver(pep_rate,  pep_on_day,  pep_off_day,  pep_active),
      prep_effective_day  = make_driver(prep_eff,  prep_on_day, prep_off_day, prep_active)
    )

    cfg_key_map <- list(
      iso_rate_day = "iso_rate_daily",
      quar_rate_day = "quar_rate_daily",
      leave_quar_day = "leave_quar_daily",
      leave_quar_prep_day = "leave_quar_prep_daily",
      prep_start_rate_day = "prep_start_daily",
      prep_stop_rate_day  = "prep_stop_daily",
      pep_rate_day = "pep_rate_daily",
      prep_effective_day  = "prep_eff_daily"
    )

    for (nm in needed_drivers)
      cfg[[ cfg_key_map[[nm]] ]] <- driver_map[[nm]]

    cfg
  }
}


################################################################################
# build_model
#
# Main entry point. Validates, generates, compiles, and returns a model object.
#
# Returns a list with:
#   $setup_sim - pre-configured setup function for this model
#   $compartments — the compartment list you passed in
#   $reactions — the reaction list you passed in
#   $cpp_file — path to the generated C++ file (for inspection)
#   $cpp_source — the generated C++ source as a string
################################################################################
#' Make a custom compartmental model from a compartment and transition rates
#'
#' @export
build_model <- function(
    compartments,
    infectious,
    in_denominator,
    reactions,
    template_path = system.file("reservoirsim_template.cpp", package = "Reservoir"),
    verbose       = TRUE
) {

  # 1. Validate
  if (verbose) message("Validating model specification...")
  validate_spec(compartments, infectious, in_denominator, reactions)

  # 2. Generate C++ source
  if (verbose) message("Generating C++ source...")
  cpp_source <- generate_cpp(
    compartments, infectious, in_denominator, reactions, template_path
  )

  # 3. Write to temp file
  cpp_file <- tempfile(pattern = "reservoirsim_", fileext = ".cpp")
  writeLines(cpp_source, cpp_file)
  if (verbose) message("Written to: ", cpp_file)

  # 4. Compile and load
  if (verbose) {
    message("Compiling model: [",
            paste(compartments[compartments != "cumul"], collapse=" -> "),
            "]")
  }
  Rcpp::sourceCpp(cpp_file)
  if (verbose) message("Done. simulate_structured_atl_cpp() is ready.")

  # 5. Return model object
  list(
    setup_sim     = make_setup_sim(compartments, reactions),
    compartments  = compartments,
    infectious    = infectious,
    in_denominator = in_denominator,
    reactions     = reactions,
    cpp_file      = cpp_file,
    cpp_source    = cpp_source
  )
}


#
# Convenience wrappers for common model structures
#

# Simple SEIR - no asymptomatic, no interventions
# Convenience constructor for a basic SEIR model
#
#' @export
model_SEIR <- function(template_path = system.file("reservoirsim_template.cpp", package = "Reservoir")) {
  build_model(
    compartments   = c("S", "L", "Is", "R", "cumul"),
    infectious     = list(list(comp="Is", weight=1.0)),
    in_denominator = c("S", "L", "Is", "R"),
    reactions = list(
      list(from="S",  to=c("L","cumul"), type="infection"),
      list(from="L",  to="Is", rate="sigma"),
      list(from="Is", to="R",  rate="gamma"),
      list(from="R",  to="S",  rate="omega"),
      list(from=NULL, to=NULL, type="birth"),
      list(from=NULL, to=NULL, type="death")
    ),
    template_path = template_path
  )
}

# Full SEIAR with isolation - matches the current human2human.cpp structure
#' Convenience constructor for an SEIAR model with isolation
#'
#' @export
model_SEIAR <- function(template_path = system.file("reservoirsim_template.cpp", package = "Reservoir")) {
  build_model(
    compartments   = c("S", "L", "Is", "Ia", "Iso", "R", "cumul"),
    infectious     = list(list(comp="Is", weight=1.0),
                          list(comp="Ia", weight="asymp_trans")),
    in_denominator = c("S", "L", "Is", "Ia", "Iso", "R"),
    reactions = list(
      list(from="S",   to=c("L","cumul"), type="infection"),
      list(from="L",   to="Is",  rate="sigma*(1-p_asym)"),
      list(from="L",   to="Ia",  rate="sigma*p_asym"),
      list(from="Is",  to="Iso", rate="iso_rate_day"),
      list(from="Is",  to="R",   rate="gamma"),
      list(from="Ia",  to="Is",  rate="kappa"),
      list(from="Ia",  to="R",   rate="gamma"),
      list(from="Iso", to="R",   rate="gamma"),
      list(from="R",   to="S",   rate="omega"),
      list(from=NULL,  to=NULL,  type="birth"),
      list(from=NULL,  to=NULL,  type="death")
    ),
    template_path = template_path
  )
}

# SEIAR with hospitalization
#' Convenience constructor for an SEIAR model with explicit hospitalization
#'
#' @export
model_SEIARH <- function(template_path = system.file("reservoirsim_template.cpp", package = "Reservoir")) {
  build_model(
    compartments   = c("S", "L", "Is", "Ia", "H", "R", "cumul"),
    infectious     = list(list(comp="Is", weight=1.0),
                          list(comp="Ia", weight="asymp_trans")),
    in_denominator = c("S", "L", "Is", "Ia", "H", "R"),
    reactions = list(
      list(from="S",  to=c("L","cumul"), type="infection"),
      list(from="L",  to="Is", rate="sigma*(1-p_asym)"),
      list(from="L",  to="Ia", rate="sigma*p_asym"),
      list(from="Is", to="H",  rate="hosp_rate"),
      list(from="Is", to="R",  rate="gamma"),
      list(from="Ia", to="Is", rate="kappa"),
      list(from="Ia", to="R",  rate="gamma"),
      list(from="H",  to="R",  rate="discharge_rate"),
      list(from="R",  to="S",  rate="omega"),
      list(from=NULL, to=NULL, type="birth"),
      list(from=NULL, to=NULL, type="death")
    ),
    template_path = template_path
  )
}
