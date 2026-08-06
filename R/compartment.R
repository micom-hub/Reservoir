################################################################################
# compartment - read one or more compartments out of a run
#
# The state matrix stores all 13 compartments for all P*A strata side by side,
# which makes manual indexing error-prone. This is the friendly way in.
################################################################################
#' Extract compartment time series from a run
#'
#' @param run A run from `generate_dataset()`.
#' @param name Compartment name(s): "S", "Sprep", "Spost", "Qs", "L", "QL",
#'   "Is", "Ia", "Iso", "R", "R2", "R3", "cumul". Vector-borne runs also have
#'   "Sv", "Ev", "Iv", "cumul_v"; waterborne runs have "W".
#' @param aggregate TRUE sums across all (patch, age) strata and returns one
#'   value per day. FALSE returns one column per stratum, named like "p0_a1".
#'
#' @return A numeric vector (one name, aggregated), or a matrix.
#'
#' @examples
#' \dontrun{
#' plot(compartment(run, "Is"), type = "l")            # infectious over time
#' matplot(compartment(run, "L", aggregate = FALSE))   # latent per stratum
#' head(compartment(run, c("S", "Is", "R")))           # several at once
#' }
#' @export
compartment <- function(run, name, aggregate = TRUE) {

  if (is.null(run$state_over_time))
    stop("This run has no compartment data. It was probably generated with ",
         "scope = \"cases_only\"; use scope = \"full\" or \"signals\".",
         call. = FALSE)

  cn  <- run$cfg$comp_names
  sot <- run$state_over_time
  P   <- run$cfg$P; A <- run$cfg$A
  blk <- P * A

  bad <- setdiff(name, cn)
  if (length(bad))
    stop("No compartment called ", paste0("\"", bad, "\"", collapse = ", "),
         ". This run has: ", paste(cn, collapse = ", "), call. = FALSE)

  one <- function(nm) {
    i <- match(nm, cn)
    # human blocks are P*A wide; appended scalars (Sv, Ev, Iv, W, ...) are 1 wide
    n_human <- run$cfg$n_human_blocks %||% 13L
    if (i <= n_human) {
      cols <- ((i - 1L) * blk + 1L):(i * blk)
    } else {
      cols <- n_human * blk + (i - n_human)
    }
    if (max(cols) > ncol(sot))
      stop("Compartment \"", nm, "\" is not present in this run's state matrix.",
           call. = FALSE)
    m <- sot[, cols, drop = FALSE]
    if (aggregate || ncol(m) == 1L) return(rowSums(m))
    colnames(m) <- paste0("p", rep(0:(P - 1), each = A),
                          "_a", rep(0:(A - 1), times = P))
    m
  }

  if (length(name) == 1L) return(one(name))

  if (!aggregate)
    stop("aggregate = FALSE works with one compartment at a time.", call. = FALSE)
  out <- vapply(name, one, numeric(nrow(sot)))
  colnames(out) <- name
  out
}
