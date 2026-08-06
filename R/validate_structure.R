################################################################################
# validate_structure.R
#
# Run it on the first few runs after any engine change.
################################################################################

.OCC_TOL <- 1e-9


################################################################################
# .block_max — peak occupancy of a named compartment across the whole run
################################################################################
.block_max <- function(state_over_time, cfg, comp) {
  idx <- match(comp, cfg$comp_names)
  if (is.na(idx)) return(NA_real_)
  block <- cfg$P * cfg$A
  cols  <- ((idx - 1L) * block + 1L):(idx * block)
  if (max(cols) > ncol(state_over_time)) return(NA_real_)
  max(state_over_time[, cols, drop = FALSE])
}


################################################################################
# validate_structure
#
# @param run   a single run object from generate_dataset()
# @param strict  TRUE stops on the first failure; FALSE returns a report
#
# Returns a data.frame of checks with columns: check, expected, observed, pass
################################################################################
#' Verify a run's realized dynamics match its declared structure
#'
#' @export
validate_structure <- function(run, strict = FALSE, verbose = TRUE) {

  cfg <- run$cfg
  s   <- run$params$structure
  sot <- run$state_over_time

  if (is.null(sot))
    stop("run has no `state_over_time`; re-run generate_dataset() with ",
         "scope = 'full' or 'signals'")

  checks <- list()
  add <- function(name, expected, observed, pass) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, expected = expected,
      observed = observed, pass = pass, stringsAsFactors = FALSE)
  }

  mx <- function(comp) .block_max(sot, cfg, comp)

  # Structural: latent bypass
  if (!s$has_latent) {
    add("L empty (no_latent)",  "0", format(mx("L"),  digits = 4), mx("L")  <= .OCC_TOL)
    add("QL empty (no_latent)", "0", format(mx("QL"), digits = 4), mx("QL") <= .OCC_TOL)
  } else {
    add("L populated (latent)", "> 0", format(mx("L"), digits = 4), mx("L") > .OCC_TOL)
  }

  # Structural: SIS
  if (s$immunity == "none") {
    add("R empty (SIS)", "0", format(mx("R"), digits = 4), mx("R") <= .OCC_TOL)
  }

  # Parametric: asymptomatic
  if (!s$has_asymptomatic)
    add("Ia empty (no asymp)", "0", format(mx("Ia"), digits = 4), mx("Ia") <= .OCC_TOL)

  # Parametric: quarantine
  if (!s$has_quarantine) {
    add("Qs empty (no quar)", "0", format(mx("Qs"), digits = 4), mx("Qs") <= .OCC_TOL)
    add("QL empty (no quar)", "0", format(mx("QL"), digits = 4), mx("QL") <= .OCC_TOL)
  }

  # Parametric: isolation
  # NOTE: Iso can also be fed by QL -> Iso (a quarantined latent who develops
  # symptoms), which is NOT gated by iso_rate. So "no isolation" only implies
  # an empty Iso when quarantine is also off, or there is no latent stage.
  if (!s$has_isolation) {
    if (!s$has_quarantine || !s$has_latent) {
      add("Iso empty (no iso)", "0", format(mx("Iso"), digits = 4), mx("Iso") <= .OCC_TOL)
    } else {
      add("Iso via QL only (no iso, quar on)", "informational",
          format(mx("Iso"), digits = 4), TRUE)
    }
  }

  # Parametric: PEP
  if (!s$has_pep)
    add("Spost empty (no PEP)", "0", format(mx("Spost"), digits = 4),
        mx("Spost") <= .OCC_TOL)

  # Parametric: vaccination / PrEP
  # This is the check that would have caught the v0.1.0 wiring bug.
  if (!s$has_vaccination) {
    add("Sprep empty (no vacc)", "0", format(mx("Sprep"), digits = 4),
        mx("Sprep") <= .OCC_TOL)
  } else if (isTRUE(run$params$has_intervention) &&
             isTRUE(run$params$prep_start_rate > 0)) {
    add("Sprep populated (vacc on)", "> 0", format(mx("Sprep"), digits = 4),
        mx("Sprep") > .OCC_TOL)
  }

  # Parametric: demography
  if (!s$has_demography) {
    add("birth_rate = 0", "0", format(cfg$params$birth_rate), cfg$params$birth_rate == 0)
    add("death_rate = 0", "0", format(cfg$params$death_rate), cfg$params$death_rate == 0)
  }

  # Parametric: recovery
  if (!s$has_recovery)
    add("gamma = 0 (SI)", "0", format(cfg$params$gamma), cfg$params$gamma == 0)

  # Immunity fate
  if (s$immunity == "permanent")
    add("omega = 0 (permanent)", "0", format(cfg$params$omega), cfg$params$omega == 0)

  # ENGINE INVARIANT: cumulative infection conservation
  # daily_new_cases is computed as a cumul-delta inside the engine. If the
  # reaction/propensity ordering were misaligned, infection events would stop
  # matching the cumul counter and this check would fail.
  cum <- run$cumul
  if (!is.null(cum) && length(cum) >= 2) {
    total_daily <- sum(run$daily_true)
    total_cumul <- cum[length(cum)] - cum[1]
    rel <- abs(total_daily - total_cumul) / max(total_cumul, 1)
    add("cumul conservation", "sum(daily) == cumul[T]-cumul[0]",
        sprintf("%.0f vs %.0f (rel %.2e)", total_daily, total_cumul, rel),
        rel < 1e-6)
  }

  # Reactive controller sanity
  if (isTRUE(cfg$reactive_mode == 1L) && !is.null(run$intervention_active)) {
    ia <- run$intervention_active
    add("reactive state is binary", "all 0/1",
        sprintf("%d distinct", length(unique(ia))),
        all(ia %in% c(0, 1)))
  }

  out <- do.call(rbind, checks)

  if (verbose) {
    skel <- s$skeleton %||% "<unknown>"
    cat(sprintf("\n--- validate_structure: %s (%s) ---\n", skel, run$mechanism))
    for (i in seq_len(nrow(out)))
      cat(sprintf("  [%s] %-38s expected %-14s got %s\n",
                  if (out$pass[i]) "OK  " else "FAIL",
                  out$check[i], out$expected[i], out$observed[i]))
    n_fail <- sum(!out$pass)
    cat(sprintf("  %d/%d checks passed\n", nrow(out) - n_fail, nrow(out)))
  }

  if (strict && any(!out$pass))
    stop("validate_structure failed: ",
         paste(out$check[!out$pass], collapse = "; "))

  invisible(out)
}


################################################################################
# validate_dataset: run validate_structure over a whole dataset
################################################################################
#' Validate every run in a generated dataset
#'
#' @export
validate_dataset <- function(dataset, max_runs = 10L, strict = FALSE) {
  all_out <- list()
  for (mech in names(dataset)) {
    runs <- dataset[[mech]]
    n <- min(length(runs), max_runs)
    if (n == 0) next
    for (i in seq_len(n)) {
      res <- validate_structure(runs[[i]], strict = strict, verbose = FALSE)
      res$mechanism <- mech
      res$run_id <- i
      all_out[[length(all_out) + 1L]] <- res
    }
  }
  if (!length(all_out)) {
    message("No runs to validate.")
    return(invisible(NULL))
  }
  out <- do.call(rbind, all_out)

  agg <- aggregate(pass ~ check, data = out,
                   FUN = function(p) sprintf("%d/%d", sum(p), length(p)))
  names(agg)[2] <- "passed"
  fails <- out[!out$pass, c("mechanism", "run_id", "check", "expected", "observed")]

  cat("\n=== validate_dataset ===\n")
  print(agg, row.names = FALSE)
  if (nrow(fails)) {
    cat("\nFailures:\n")
    print(fails, row.names = FALSE)
  } else {
    cat("\nAll checks passed.\n")
  }
  invisible(list(all = out, failures = fails))
}
