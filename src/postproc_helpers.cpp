//
// postproc_helpers.cpp
//
// Two helpers for postprocessing.R
//

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;


//
// apply_reporting_delays_cpp
//
// @param cases            daily case counts (length n)
// @param progress         0-1 outbreak progress per day (length n); drives the
//                         logistic shrink of the delay distribution
// @param init_max_delay   maximum reporting delay early in the outbreak (days)
// @param final_max_delay  maximum reporting delay late in the outbreak (days)
// @param gamma_shape      shape of the within-window delay kernel
// @param gamma_rate    rate  of the within-window delay kernel
// @param steepness  logistic steepness for the shrink
// @param midpoint   logistic midpoint (on the progress scale)
//
//
// [[Rcpp::export]]
NumericVector apply_reporting_delays_cpp(NumericVector cases,
                                         NumericVector progress,
                                         int init_max_delay,
                                         int final_max_delay,
                                         double gamma_shape,
                                         double gamma_rate,
                                         double steepness,
                                         double midpoint) {
  int n = cases.size();
  NumericVector out(n, 0.0);
  if (n == 0) return out;

  if (init_max_delay  < 0) init_max_delay  = 0;
  if (final_max_delay < 0) final_max_delay = 0;
  if (gamma_shape <= 0.0) gamma_shape = 1.0;
  if (gamma_rate  <= 0.0) gamma_rate  = 1.0;

  int hard_max = std::max(init_max_delay, final_max_delay);

  // Pre-compute the unnormalized Gamma kernel once
  std::vector<double> kern(hard_max + 1, 0.0);
  for (int d = 0; d <= hard_max; ++d) {
    // density at the midpoint of day d
    double x = static_cast<double>(d) + 0.5;
    kern[d] = R::dgamma(x, gamma_shape, 1.0 / gamma_rate, 0);
  }

  for (int t = 0; t < n; ++t) {
    double c = cases[t];
    if (c <= 0.0) continue;

    // Logistic shrink of the maximum delay as the outbreak progresses.
    double p = (t < progress.size()) ? progress[t] : 1.0;
    if (!R_finite(p)) p = 0.0;
    double w = 1.0 / (1.0 + std::exp(-steepness * (p - midpoint)));
    int max_delay = static_cast<int>(std::lround(
      init_max_delay + w * (final_max_delay - init_max_delay)));
    if (max_delay < 0)        max_delay = 0;
    if (max_delay > hard_max) max_delay = hard_max;

    // Normalize the kernel over the active window.
    double ksum = 0.0;
    for (int d = 0; d <= max_delay; ++d) ksum += kern[d];
    if (ksum <= 0.0) { out[t] += c; continue; }

    for (int d = 0; d <= max_delay; ++d) {
      double share = c * kern[d] / ksum;
      int tgt = t + d;
      if (tgt >= n) tgt = n - 1;   // retain tail mass rather than discard it
      out[tgt] += share;
    }
  }

  return out;
}


//
// hospital_occupancy_cpp
//
// @param daily_hosp  daily admissions (length n)
// @param avg_los     mean length of stay in days
//
// [[Rcpp::export]]
IntegerVector hospital_occupancy_cpp(IntegerVector daily_hosp, int avg_los) {
  int n = daily_hosp.size();
  IntegerVector occ(n);   // Rcpp zero-initialises
  if (n == 0) return occ;
  if (avg_los < 1) avg_los = 1;

  double shape = 2.0;
  double scale = static_cast<double>(avg_los) / shape;

  // Above this many admissions in a day, sample a mean LOS for the whole batch instead of per-individual
  const int BATCH_THRESHOLD = 500;

  for (int t = 0; t < n; ++t) {
    int adm = daily_hosp[t];
    if (adm <= 0) continue;

    if (adm <= BATCH_THRESHOLD) {
      for (int k = 0; k < adm; ++k) {
        int los = static_cast<int>(std::lround(R::rgamma(shape, scale)));
        if (los < 1) los = 1;
        int end = std::min(t + los, n);
        for (int d = t; d < end; ++d) occ[d] += 1;
      }
    } else {
      // Batched: distribute admissions across Gamma LOS.
      int max_los = std::max(1, static_cast<int>(std::lround(avg_los * 5.0)));
      std::vector<double> w(max_los + 1, 0.0);
      double wsum = 0.0;
      for (int L = 1; L <= max_los; ++L) {
        w[L] = R::dgamma(static_cast<double>(L), shape, scale, 0);
        wsum += w[L];
      }
      if (wsum <= 0.0) { w[avg_los] = 1.0; wsum = 1.0; }

      for (int L = 1; L <= max_los; ++L) {
        double share = adm * w[L] / wsum;
        int n_pat = static_cast<int>(std::lround(share));
        if (n_pat <= 0) continue;
        int end = std::min(t + L, n);
        for (int d = t; d < end; ++d) occ[d] += n_pat;
      }
    }
  }

  return occ;
}
