###############################################################################
# generate_slots.R
###############################################################################

###############################################################################
# Helper: capitalize first letter of a string
# "cumul" -> "Cumul",  "Is" -> "Is",  "S" -> "S"
# Used so all C++ variable names follow UpperCamelCase after the prefix:
#   offCumul, iCumul  (not offcumul, icumul)
###############################################################################
cap1 <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}


###############################################################################
# SLOT 1 - N_BLOCKS
###############################################################################
gen_n_blocks <- function(compartments) {
  as.character(length(compartments))
}


###############################################################################
# SLOT 2 - OFFSET_PARAMS
# int offS, int offL, int offIs, ... int offCumul
###############################################################################
gen_offset_params <- function(compartments) {
  params <- paste0("    int off", cap1(compartments))
  paste(params, collapse = ",\n")
}


###############################################################################
# SLOT 3 - LOCAL_INDEX_VARS  (inside build_reactions loop)
# int iS = idx_of(p, a, offS);
# int iCumul = idx_of(p, a, offCumul);
###############################################################################
gen_local_index_vars <- function(compartments) {
  max_len <- max(nchar(cap1(compartments)))
  lines   <- vapply(compartments, function(c) {
    cc  <- cap1(c)
    pad <- paste(rep(" ", max_len - nchar(cc)), collapse = "")
    sprintf("int i%s%s = idx_of(p, a, off%s);", cc, pad, cc)
  }, character(1))
  paste(lines, collapse = "\n")
}


###############################################################################
# SLOT 4 — REACTION_LINES  (inside build_reactions loop)
###############################################################################
gen_reaction_lines <- function(reactions, compartments) {
  mortal <- compartments[compartments != "cumul"]

  lines <- vapply(reactions, function(rxn) {
    type <- if (!is.null(rxn$type)) rxn$type else "flow"

    if (type == "infection") {
      targets <- rxn$to
      all_c   <- c(rxn$from, targets)
      deltas  <- c(-1L, rep(+1L, length(targets)))
      idx_str <- paste0("i", cap1(all_c), collapse = ", ")
      del_str <- paste(deltas, collapse = ", ")
      return(sprintf("add({%s}, {%s});  // infection: %s -> %s",
                     idx_str, del_str, rxn$from,
                     paste(targets, collapse = "+")))
    }

    if (type == "birth") {
      first <- cap1(mortal[1])
      return(sprintf("add({i%s}, {+1});  // birth", first))
    }

    if (type == "death") {
      dl <- vapply(mortal, function(c)
        sprintf("add({i%s}, {-1});", cap1(c)), character(1))
      return(paste(dl, collapse = "\n      "))
    }

    # Standard flow
    sprintf("add({i%s, i%s}, {-1, +1});  // %s -> %s  rate: %s",
            cap1(rxn$from), cap1(rxn$to),
            rxn$from, rxn$to,
            if (!is.null(rxn$rate)) rxn$rate else "?")

  }, character(1))

  paste(lines, collapse = "\n      ")
}


###############################################################################
# SLOT 5 — LAMBDA_VAR_DECLS
# double S = get(q, b, offS);  - only for denominator + infectious compartments
###############################################################################
gen_lambda_var_decls <- function(in_denominator, infectious) {
  inf_comps <- vapply(infectious, `[[`, character(1), "comp")
  needed    <- unique(c(in_denominator, inf_comps))
  max_len   <- max(nchar(needed))

  lines <- vapply(needed, function(c) {
    pad <- paste(rep(" ", max_len - nchar(c)), collapse = "")
    sprintf("double %s%s = get(q, b, off%s);", c, pad, cap1(c))
  }, character(1))

  paste(lines, collapse = "\n          ")
}


###############################################################################
# SLOT 6 — N_SUM
# "S + L + Is + Ia + Iso + R"
###############################################################################
gen_n_sum <- function(in_denominator) {
  paste(in_denominator, collapse = " + ")
}


################################################################################
# SLOT 7 — IEFF
# "Is + asymp_trans * Ia"
###############################################################################
gen_ieff <- function(infectious) {
  if (length(infectious) == 0) return("0.0")
  terms <- vapply(infectious, function(i) {
    w <- i$weight
    if ((is.numeric(w) && w == 1.0) || identical(w, "1") || identical(w, "1.0"))
      i$comp
    else
      paste0(w, " * ", i$comp)
  }, character(1))
  paste(terms, collapse = " + ")
}


###############################################################################
# SLOT 8 — PROP_VAR_DECLS  (inside compute_propensities loop)
# double S = get(p, aa, offS);  — all non-cumul compartments
###############################################################################
gen_prop_var_decls <- function(compartments) {
  mortal  <- compartments[compartments != "cumul"]
  max_len <- max(nchar(mortal))
  lines   <- vapply(mortal, function(c) {
    pad <- paste(rep(" ", max_len - nchar(c)), collapse = "")
    sprintf("double %s%s = get(p, aa, off%s);", c, pad, cap1(c))
  }, character(1))
  paste(lines, collapse = "\n      ")
}


###############################################################################
# SLOT 9 — PROPENSITY_LINES
# Same iteration order as gen_reaction_lines — this is the ordering guarantee.
###############################################################################
gen_propensity_lines <- function(reactions, compartments) {
  mortal <- compartments[compartments != "cumul"]

  lines <- vapply(reactions, function(rxn) {
    type <- if (!is.null(rxn$type)) rxn$type else "flow"

    if (type == "infection") {
      modifier <- rxn$modifier
      if (is.null(modifier))
        return(sprintf("a[j++] = lam * %s;  // infection from %s",
                       rxn$from, rxn$from))
      else
        return(sprintf("a[j++] = lam * (%s) * %s;  // infection from %s",
                       modifier, rxn$from, rxn$from))
    }

    if (type == "birth")
      return("a[j++] = birth_rate * N;  // birth")

    if (type == "death") {
      dl <- vapply(mortal, function(c)
        sprintf("a[j++] = death_rate * %s;", c), character(1))
      return(paste(dl, collapse = "\n      "))
    }

    rate <- if (!is.null(rxn$rate)) rxn$rate else
      stop("Reaction from='", rxn$from, "' to='", rxn$to, "' has no rate")

    sprintf("a[j++] = %s * %s;  // %s -> %s",
            rate, rxn$from, rxn$from, rxn$to)

  }, character(1))

  paste(lines, collapse = "\n      ")
}


###############################################################################
# SLOT 10 — OFFSET_ASSIGNMENTS  (in main function)
# int offS     = 0L * block;
# int offCumul = 6L * block;
###############################################################################
gen_offset_assignments <- function(compartments) {
  max_len <- max(nchar(cap1(compartments)))
  lines   <- vapply(seq_along(compartments), function(i) {
    cc  <- cap1(compartments[i])
    pad <- paste(rep(" ", max_len - nchar(cc)), collapse = "")
    sprintf("int off%s%s = %dL * block;", cc, pad, i - 1L)
  }, character(1))
  paste(lines, collapse = "\n  ")
}


###############################################################################
# SLOT 11 — OFFSET_ARGS  (call sites)
# offS, offL, offIs, ... offCumul
###############################################################################
gen_offset_args <- function(compartments) {
  args <- paste0("    off", cap1(compartments))
  paste(args, collapse = ",\n")
}


###############################################################################
# SLOT 12 — CUMUL_OFFSET_IDX
###############################################################################
gen_cumul_offset_idx <- function(compartments) {
  idx <- which(tolower(compartments) == "cumul") - 1L
  if (length(idx) == 0) stop("compartments must include 'cumul'")
  sprintf("%dL * block", idx)
}


###############################################################################
# Driver detection -  scans rate strings for known driver tokens
###############################################################################
.DRIVER_MAP <- list(
  iso_rate_day        = list(cfg = "iso_rate_daily",        var = "iso_rate_day"),
  quar_rate_day       = list(cfg = "quar_rate_daily",       var = "quar_rate_day"),
  leave_quar_day      = list(cfg = "leave_quar_daily",      var = "leave_quar_day"),
  leave_quar_prep_day = list(cfg = "leave_quar_prep_daily", var = "leave_quar_prep_day"),
  prep_start_rate_day = list(cfg = "prep_start_daily",      var = "prep_start_rate_day"),
  prep_stop_rate_day  = list(cfg = "prep_stop_daily",       var = "prep_stop_rate_day"),
  pep_rate_day        = list(cfg = "pep_rate_daily",        var = "pep_rate_day"),
  prep_effective_day  = list(cfg = "prep_eff_daily",        var = "prep_effective_day")
)

detect_needed_drivers <- function(reactions) {
  # Collect all rate and modifier strings safely
  rate_strs <- c()
  for (r in reactions) {
    if (!is.null(r$rate))     rate_strs <- c(rate_strs, r$rate)
    if (!is.null(r$modifier)) rate_strs <- c(rate_strs, r$modifier)
  }
  if (length(rate_strs) == 0) return(character(0))
  all_text <- paste(rate_strs, collapse = " ")
  Filter(function(nm) grepl(nm, all_text, fixed = TRUE), names(.DRIVER_MAP))
}

# Extra biological params: tokens in rate strings not already in the fixed set
.FIXED_TOKENS <- c("sigma","gamma","omega","p_asym","asymp_trans","kappa",
                   "birth_rate","death_rate","lam",
                   "iso_rate_day","quar_rate_day","leave_quar_day",
                   "leave_quar_prep_day","prep_start_rate_day",
                   "prep_stop_rate_day","pep_rate_day","prep_effective_day")

detect_extra_params <- function(reactions) {
  rate_strs <- c()
  for (r in reactions) {
    if (!is.null(r$rate))     rate_strs <- c(rate_strs, r$rate)
    if (!is.null(r$modifier)) rate_strs <- c(rate_strs, r$modifier)
  }
  if (length(rate_strs) == 0) return(character(0))
  tokens <- unique(unlist(regmatches(
    rate_strs, gregexpr("[A-Za-z_][A-Za-z0-9_]*", rate_strs)
  )))
  setdiff(tokens, .FIXED_TOKENS)
}


# -----------------------------------------------------------------------------
# SLOT 13 — DRIVER_READS + DRIVER_DRV_LINES
# -----------------------------------------------------------------------------
gen_driver_reads <- function(reactions) {
  needed <- detect_needed_drivers(reactions)
  if (length(needed) == 0) return("  // no time-varying intervention drivers")
  lines <- vapply(needed, function(nm) {
    d <- .DRIVER_MAP[[nm]]
    sprintf('  NumericVector %s_vec = as<NumericVector>(cfg["%s"]);',
            nm, d$cfg)
  }, character(1))
  paste(lines, collapse = "\n")
}

gen_driver_drv_lines <- function(reactions) {
  needed <- detect_needed_drivers(reactions)
  if (length(needed) == 0) return("    // no intervention drivers")
  lines <- vapply(needed, function(nm)
    sprintf("    double %s = drv(%s_vec, day);", nm, nm),
    character(1))
  paste(lines, collapse = "\n")
}


###############################################################################
# SLOT 14 — PARAM_READS
###############################################################################
gen_param_reads <- function(reactions) {
  fixed <- c(
    '  double sigma       = as<double>(par["sigma"]);',
    '  double gamma       = as<double>(par["gamma"]);',
    '  double omega       = as<double>(par["omega"]);',
    '  double p_asym      = as<double>(par["p_asym"]);',
    '  double asymp_trans = as<double>(par["asymp_trans"]);',
    '  double kappa       = as<double>(par["kappa"]);',
    '  double birth_rate  = as<double>(par["birth_rate"]);',
    '  double death_rate  = as<double>(par["death_rate"]);'
  )
  extra <- detect_extra_params(reactions)
  extra_lines <- vapply(extra, function(tok)
    sprintf('  double %s = as<double>(par["%s"]);', tok, tok),
    character(1))
  paste(c(fixed, extra_lines), collapse = "\n")
}


###############################################################################
# SLOT 15 — RATE_PARAMS + RATE_ARGS
###############################################################################
gen_rate_params <- function(reactions) {
  needed <- detect_needed_drivers(reactions)
  extra  <- detect_extra_params(reactions)

  fixed_params <- c("sigma","gamma","omega","p_asym","asymp_trans",
                    "kappa","birth_rate","death_rate")
  all_params <- c(fixed_params, extra, needed)
  paste(paste0("    double ", all_params), collapse = ",\n")
}

gen_rate_args <- function(reactions) {
  needed <- detect_needed_drivers(reactions)
  extra  <- detect_extra_params(reactions)

  fixed <- c("sigma","gamma","omega","p_asym","asymp_trans",
             "kappa","birth_rate","death_rate")
  all_args <- c(fixed, extra, needed)
  paste(paste0("    ", all_args), collapse = ",\n")
}


###############################################################################
# fill_template
###############################################################################
fill_template <- function(template_path, slots) {
  txt <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
  for (nm in names(slots)) {
    txt <- gsub(paste0("%%", nm, "%%"), slots[[nm]], txt, fixed = TRUE)
  }
  remaining <- regmatches(txt, gregexpr("%%[A-Z_]+%%", txt))[[1]]
  if (length(remaining) > 0)
    warning("Unfilled placeholders: ", paste(remaining, collapse = ", "))
  txt
}


###############################################################################
# test_generators — run after source() to verify output
###############################################################################
test_generators <- function() {

  comps <- c("S", "L", "Is", "Ia", "Iso", "R", "cumul")
  infec <- list(list(comp="Is", weight=1.0),
                list(comp="Ia", weight="asymp_trans"))
  denom <- c("S", "L", "Is", "Ia", "Iso", "R")

  rxns <- list(
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
  )

  cat("--- N_BLOCKS ---\n");          cat(gen_n_blocks(comps), "\n\n")
  cat("--- OFFSET_PARAMS ---\n");     cat(gen_offset_params(comps), "\n\n")
  cat("--- LOCAL_INDEX_VARS ---\n");  cat(gen_local_index_vars(comps), "\n\n")
  cat("--- REACTION_LINES ---\n");    cat(gen_reaction_lines(rxns, comps), "\n\n")
  cat("--- LAMBDA_VAR_DECLS ---\n");  cat(gen_lambda_var_decls(denom, infec), "\n\n")
  cat("--- N_SUM ---\n");             cat(gen_n_sum(denom), "\n\n")
  cat("--- IEFF ---\n");              cat(gen_ieff(infec), "\n\n")
  cat("--- PROP_VAR_DECLS ---\n");    cat(gen_prop_var_decls(comps), "\n\n")
  cat("--- PROPENSITY_LINES ---\n");  cat(gen_propensity_lines(rxns, comps), "\n\n")
  cat("--- OFFSET_ASSIGNMENTS ---\n");cat(gen_offset_assignments(comps), "\n\n")
  cat("--- OFFSET_ARGS ---\n");       cat(gen_offset_args(comps), "\n\n")
  cat("--- CUMUL_OFFSET_IDX ---\n");  cat(gen_cumul_offset_idx(comps), "\n\n")
  cat("--- DRIVER_READS ---\n");      cat(gen_driver_reads(rxns), "\n\n")
  cat("--- DRIVER_DRV_LINES ---\n");  cat(gen_driver_drv_lines(rxns), "\n\n")
  cat("--- PARAM_READS ---\n");       cat(gen_param_reads(rxns), "\n\n")
  cat("--- RATE_PARAMS ---\n");       cat(gen_rate_params(rxns), "\n\n")
  cat("--- RATE_ARGS ---\n");         cat(gen_rate_args(rxns), "\n\n")

  invisible(TRUE)
}

# source("generate_slots.R")
# test_generators()
