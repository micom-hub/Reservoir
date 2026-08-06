################################################################################
# aaa_constants.R
#
# Package-wide constants and infix helpers.
#
# Named with an `aaa_` purposefully: R sources R/*.R in collation
# (alphabetical) order at build time, so any symbol referenced by an upstream
# funciton in another file must be defined in a file that runs earlier.
# `.HUMAN_COMPS` is used further upstream by extract_compartment_csvs.R, which
# sorts before setup_sim.R. Do not rename this file!!!!!!
################################################################################

# This is the human compartment block order and it needs to be
# identical across all three transmission mechanism simulators


.HUMAN_COMPS <- c("S", "Sprep", "Spost", "Qs",
                  "L", "QL", "Is", "Ia", "Iso",
                  "R", "R2", "R3", "cumul")

.N_HUMAN_BLOCKS <- 13

# NULL-coalescing helper used throughout the package.
`%||%` <- function(a,b) if (is.null(a)) b else a


################################################################################
# .sample1 — sample one element, safely
#
# R's sample() has a long-standing trap: sample(5, 1) draws from 1:5 rather
# than returning 5. That bites whenever a range is narrowed to a single value.
################################################################################-
.sample1 <- function(x, prob = NULL) {
  if (length(x) == 1L) return(x)
  x[sample.int(length(x), 1L, prob = prob)]
}


################################################################################
# .as_dist: this function accepts a scalar, a range, etc. for a distribution
#
#   6            pin the value  (e.g. infectious_days = 6)
#   c(3, 10)     a range; interpreted as a central ~90% interval
#   list(...)    full control, passed through unchanged
#
# `type` selects the field names the sampler expects:
#   "clipped"  meanlog, sdlog, lo, hi      (durations)
#   "shifted"  meanlog, sdlog, shift, floor (R0 priors)
################################################################################
.as_dist <- function(x, nm, type = c("clipped", "shifted")) {
  type <- match.arg(type)
  if (is.null(x)) return(NULL)
  if (is.list(x)) return(x)

  if (!is.numeric(x) || length(x) > 2L || any(!is.finite(x)) || any(x <= 0))
    stop("`", nm, "` must be a positive number, a length-2 range c(lo, hi), ",
         "or a full list(meanlog=, sdlog=, ...)")

  if (length(x) == 1L) {
    return(if (type == "shifted")
      list(meanlog = log(x), sdlog = 1e-6, shift = 0, floor = x * 0.999)
    else
      list(meanlog = log(x), sdlog = 1e-6, lo = x * 0.999, hi = x * 1.001))
  }

  lo <- min(x); hi <- max(x)
  if (lo == hi) return(.as_dist(lo, nm, type))
  meanlog <- log(sqrt(lo * hi))
  sdlog   <- log(hi / lo) / (2 * 1.6449)
  if (type == "shifted")
    list(meanlog = meanlog, sdlog = sdlog, shift = 0, floor = lo)
  else
    list(meanlog = meanlog, sdlog = sdlog, lo = lo, hi = hi)
}
