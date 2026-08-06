################################################################################
# extract_compartment_csvs.R
#
# Exports saved runs to CSV
#
################################################################################

.MECH_COMPS <- list(
  base       = .HUMAN_COMPS,
  vector     = c(.HUMAN_COMPS, "Sv", "Ev", "Iv", "cumul_v"),
  waterborne = c(.HUMAN_COMPS, "W")
)


.guess_mechanism <- function(path) {
  if (grepl("vector", path, ignore.case = TRUE)) return("vector")
  if (grepl("water",  path, ignore.case = TRUE)) return("waterborne")
  "base"
}


################################################################################
# .run_layout
################################################################################
.run_layout <- function(run, fallback_mech = "base") {
  cfg <- run$cfg
  comp_names <- cfg$comp_names
  if (is.null(comp_names)) {
    mech <- run$mechanism %||% fallback_mech
    comp_names <- .MECH_COMPS[[mech]] %||% .HUMAN_COMPS
  }
  n_human <- cfg$n_human_blocks %||% .N_HUMAN_BLOCKS
  list(comp_names = comp_names, n_human = n_human,
       P = cfg$P, A = cfg$A, tf = cfg$tf_days, block = cfg$P * cfg$A)
}


################################################################################
# extract_compartment_csvs
################################################################################
#' Write one CSV per run of compartment time series
#'
#' @export
extract_compartment_csvs <- function(rds_path, out_dir,
                                     mechanism     = NULL,
                                     aggregate     = TRUE,
                                     include_noise = FALSE) {

  if (!file.exists(rds_path)) stop("File not found: ", rds_path)
  runs <- readRDS(rds_path)
  if (!length(runs)) { message("No runs in ", rds_path); return(invisible(NULL)) }
  if (is.null(runs[[1]]$state_over_time))
    stop("state_over_time not found. Re-run generate_dataset() with ",
         "scope = 'full' or 'signals'.")

  if (is.null(mechanism)) mechanism <- .guess_mechanism(rds_path)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Extracting ", length(runs), " runs [", mechanism, "] -> ", out_dir,
          if (aggregate) " (aggregated)" else " (per patch/age)")

  for (i in seq_along(runs)) {
    df <- .run_to_df(runs[[i]], mechanism, aggregate, include_noise)
    write.csv(df, file.path(out_dir, sprintf("run_%03d.csv", i)), row.names = FALSE)
  }

  message("Done. ", length(runs), " csvs written to ", out_dir)
  invisible(out_dir)
}


################################################################################
# .run_to_df — convert one run to a data.frame
################################################################################
.run_to_df <- function(r, fallback_mech, aggregate, include_noise) {

  lay <- .run_layout(r, fallback_mech)
  sot <- r$state_over_time
  P <- lay$P; A <- lay$A; block <- lay$block
  M_human <- lay$n_human * block

  df <- data.frame(day = seq_len(lay$tf))

  for (ci in seq_along(lay$comp_names)) {
    nm <- lay$comp_names[ci]

    if (ci <= lay$n_human) {
      start_col <- (ci - 1L) * block + 1L
      end_col   <- ci * block

      if (end_col > ncol(sot)) {
        if (aggregate || block == 1L) {
          df[[nm]] <- NA_real_
        } else {
          for (p in 0:(P - 1)) for (a in 0:(A - 1))
            df[[sprintf("%s_p%d_a%d", nm, p, a)]] <- NA_real_
        }
      } else if (aggregate || block == 1L) {
        df[[nm]] <- if (block == 1L) sot[, start_col]
        else rowSums(sot[, start_col:end_col, drop = FALSE])
      } else {
        # layout within a block: patch varies slowly, age varies fast
        for (p in 0:(P - 1)) for (a in 0:(A - 1)) {
          col_idx <- (ci - 1L) * block + p * A + a + 1L
          df[[sprintf("%s_p%d_a%d", nm, p, a)]] <- sot[, col_idx]
        }
      }

    } else {
      # scalar slots appended after the human blocks
      col_idx <- M_human + (ci - lay$n_human)
      df[[nm]] <- if (col_idx > ncol(sot)) NA_real_ else sot[, col_idx]
    }
  }

  # True incidence and Rt are core outputs, not optional extras: every export
  # carries the case series. include_noise controls only the derived and
  # noised surveillance streams.
  df$cases_true <- r$daily_true
  df$Rt <- if (!is.null(r$daily_Rt)) r$daily_Rt[seq_len(nrow(df))] else NA_real_

  if (include_noise) {
    df$cases_limited <- if (!is.null(r$noise_limited)) r$noise_limited$daily_cases else NA

    if (!is.null(r$noise_full)) {
      nf <- r$noise_full
      df$cases_reported  <- nf$daily_cases_reported
      df$hosp_true <- nf$daily_hosp_true
      df$hosp_reported <- nf$daily_hosp_reported
      df$deaths_true <- nf$daily_death_true
      df$deaths_reported <- nf$daily_death_reported
      df$hosp_occupancy <- nf$hosp_occupancy
      df$wastewater_raw <- nf$wastewater_raw
      df$wastewater_norm <- nf$wastewater_normalized
      df$syndromic <- nf$syndromic_signal
    }
  }

  if (!is.null(r$intervention_active))
    df$intervention_active <- r$intervention_active[seq_len(nrow(df))]

  df
}


################################################################################
# extract_all_runs_single_csv
#
# Combines all runs into one long CSV with run_id, mechanism, Rt, and the
# full set of structure labels up front
#
################################################################################
#' Combine all runs into a single long CSV
#'
#' @export
extract_all_runs_single_csv <- function(rds_path, out_path,
                                        mechanism     = NULL,
                                        aggregate     = TRUE,
                                        include_noise = FALSE) {

  if (!file.exists(rds_path)) stop("File not found: ", rds_path)
  if (is.null(mechanism)) mechanism <- .guess_mechanism(rds_path)

  runs <- readRDS(rds_path)
  if (!length(runs)) { message("No runs in ", rds_path); return(invisible(NULL)) }

  all_dfs <- lapply(seq_along(runs), function(i) {
    r  <- runs[[i]]
    df <- .run_to_df(r, mechanism, aggregate, include_noise)

    df$run_id <- i
    df$mechanism <- r$mechanism %||% mechanism

    # structure labels (constant within a run, broadcast across days)
    s <- r$params$structure
    if (!is.null(s)) {
      df$skeleton <- s$skeleton
      df$has_latent <- s$has_latent
      df$has_asymptomatic <- s$has_asymptomatic
      df$has_recovery <- s$has_recovery
      df$immunity <- s$immunity
      df$n_waning_stages <- s$n_waning_stages
      df$has_isolation <- s$has_isolation
      df$has_quarantine <- s$has_quarantine
      df$has_pep <- s$has_pep
      df$has_vaccination <- s$has_vaccination
      df$has_demography <- s$has_demography
      df$has_importation <- s$has_importation
      df$has_seasonality <- isTRUE(r$params$use_seasonality)
    }
    df$is_vectorborne <- as.integer((r$mechanism %||% mechanism) == "vector")
    df$is_waterborne <- as.integer((r$mechanism %||% mechanism) == "waterborne")
    df$population <- r$cfg$N
    df
  })

  # Union columns across runs (structures may differ between runs)
  all_cols <- unique(unlist(lapply(all_dfs, names)))
  all_dfs <- lapply(all_dfs, function(d) {
    missing <- setdiff(all_cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })

  combined <- do.call(rbind, all_dfs)

  front <- c("run_id", "mechanism", "day", "cases_true", "Rt",
             "population", "skeleton",
             "has_latent", "has_asymptomatic", "has_recovery", "immunity",
             "n_waning_stages", "has_isolation", "has_quarantine", "has_pep",
             "has_vaccination", "has_demography", "has_importation",
             "has_seasonality", "is_vectorborne", "is_waterborne")
  front <- intersect(front, names(combined))
  rest <- setdiff(names(combined), front)
  combined <- combined[, c(front, rest), drop = FALSE]

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(combined, out_path, row.names = FALSE)

  message("Written: ", out_path)
  message("  ", nrow(combined), " rows  x  ", ncol(combined), " columns")
  message("  runs: ", length(runs),
          "  days per run: ", round(nrow(combined) / length(runs), 1))

  invisible(out_path)
}
