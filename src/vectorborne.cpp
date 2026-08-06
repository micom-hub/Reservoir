#include <Rcpp.h>
using namespace Rcpp;

namespace {  // file-local helpers — avoids cross-TU symbol clashes

// ---------- small helpers ----------
static inline double safe_div(double a, double b) { return (b > 0.0) ? (a / b) : 0.0; }
static inline double rexp_rate(double rate) { return (rate > 0.0) ? R::rexp(1.0 / rate) : R_PosInf; }
static inline bool is_finite(double x) { return R_finite(x); }

// Sparse reaction: indices + deltas (both same length)
struct Reaction {
  std::vector<int> idx;
  std::vector<int> delta;
};

// Apply reaction delta to state (in-place)
static inline void apply_reaction(std::vector<double>& x, const Reaction& r) {
  for (size_t k = 0; k < r.idx.size(); ++k)
    x[r.idx[k]] += static_cast<double>(r.delta[k]);
}

// Check if reaction is "critical": any consumed species below Ncritical
static inline bool is_critical(const std::vector<double>& x,
                               const Reaction& r, int Ncritical) {
  for (size_t k = 0; k < r.idx.size(); ++k)
    if (r.delta[k] < 0 && x[r.idx[k]] < static_cast<double>(Ncritical))
      return true;
    return false;
}

// ------------------------------------------------------------
// Human compartment blocks (13 total, same as human2human.cpp):
//   0  S
//   1  Sprep
//   2  Spost
//   3  Qs
//   4  L
//   5  QL
//   6  Is
//   7  Ia
//   8  Iso
//   9  R        (first waning stage, always present)
//  10  R2       (second waning stage)
//  11  R3       (third waning stage)
//  12  cumul    (cumulative infections tracker)
//
// Vector compartments (4 scalars appended after human blocks):
//   M_human+0  Sv      (susceptible vectors)
//   M_human+1  Ev      (exposed vectors)
//   M_human+2  Iv      (infectious vectors)
//   M_human+3  cumul_v (cumulative vector infections)
//
// Total state vector length = 13*P*A + 4
// ------------------------------------------------------------
static constexpr int N_HUMAN_BLOCKS = 13;

// Build reactions.
// Human reactions per (p,a): 38
//   Base flow reactions: 22  (infections, PrEP, quarantine, PEP, progression, recovery)
//   Waning:               3  (R -> R2 -> R3 -> S)
//   Birth:                1
//   Deaths:              12  (one per mortal human compartment)
// Vector reactions (global, appended once): 6
static std::vector<Reaction> build_reactions(
    int P, int A, int no_latent, int sis_mode,
    int offS,    int offSprep, int offSpost, int offQs,
    int offL,    int offQL,    int offIs,    int offIa,  int offIso,
    int offR,    int offR2,    int offR3,    int offCumul,
    int iSv, int iEv, int iIv, int iCumulV
) {
  std::vector<Reaction> rxn;
  static constexpr int RXN_PER_PA = 38;
  // Net change when the latent stage is bypassed:
  //   +2 infection channels (direct Is/Ia split)
  //   -7 L/QL reactions (L->QL, L->Spost, QL->Spost, L->Is, L->Ia,
  //                      QL->Ia, QL->Iso)  =>  net -5 per (p,a) cell
  static constexpr int RXN_LATENT_DELTA = 5;
  static constexpr int RXN_VECTOR =  6;
  rxn.reserve(static_cast<size_t>(P) * static_cast<size_t>(A) * RXN_PER_PA + RXN_VECTOR);
  
  auto idx_of = [&](int p, int a, int offset) {
    return offset + (p * A + a);
  };
  
  for (int p = 0; p < P; ++p) {
    for (int a = 0; a < A; ++a) {
      int iS     = idx_of(p, a, offS);
      int iSprep = idx_of(p, a, offSprep);
      int iSpost = idx_of(p, a, offSpost);
      int iQs    = idx_of(p, a, offQs);
      int iL     = idx_of(p, a, offL);
      int iQL    = idx_of(p, a, offQL);
      int iIs    = idx_of(p, a, offIs);
      int iIa    = idx_of(p, a, offIa);
      int iIso   = idx_of(p, a, offIso);
      int iR     = idx_of(p, a, offR);
      int iR2    = idx_of(p, a, offR2);
      int iR3    = idx_of(p, a, offR3);
      int iCum   = idx_of(p, a, offCumul);
      
      auto add = [&](std::initializer_list<int> idxs,
                     std::initializer_list<int> dels) {
        Reaction r;
        r.idx.assign(idxs.begin(), idxs.end());
        r.delta.assign(dels.begin(), dels.end());
        rxn.push_back(std::move(r));
      };
      
      // Recovery destination: R normally, S under SIS (immunity = "none")
      int iRec = sis_mode ? iS : iR;

      // ---- Infection ----
      if (!no_latent) {
        add({iS,     iL, iCum}, {-1, +1, +1});  // S     -> L
        add({iSprep, iL, iCum}, {-1, +1, +1});  // Sprep -> L (PrEP, reduced)
      } else {
        // NO-LATENT: infections deposit directly into Is / Ia, split by p_asym
        add({iS,     iIs, iCum}, {-1, +1, +1}); // S     -> Is
        add({iS,     iIa, iCum}, {-1, +1, +1}); // S     -> Ia
        add({iSprep, iIs, iCum}, {-1, +1, +1}); // Sprep -> Is
        add({iSprep, iIa, iCum}, {-1, +1, +1}); // Sprep -> Ia
      }
      
      // ---- PrEP transitions ----
      add({iS,     iSprep}, {-1, +1});
      add({iSprep, iS},     {-1, +1});
      
      // ---- Quarantine (susceptible) ----
      add({iS,     iQs},    {-1, +1});
      add({iQs,    iS},     {-1, +1});
      add({iQs,    iSprep}, {-1, +1});
      
      if (!no_latent) {
        add({iL,     iQL},    {-1, +1});        // L -> QL
      }
      
      // ---- Sprep quarantine ----
      add({iSprep, iQs},    {-1, +1});
      
      if (!no_latent) {
        add({iL,     iSpost}, {-1, +1});        // L  -> Spost (PEP)
        add({iQL,    iSpost}, {-1, +1});        // QL -> Spost (PEP)
      }
      add({iSpost, iS},     {-1, +1});          // Spost -> S
      
      if (!no_latent) {
        add({iL,     iIs},    {-1, +1});
        add({iL,     iIa},    {-1, +1});
        add({iQL,    iIa},    {-1, +1});
        add({iQL,    iIso},   {-1, +1});
      }
      
      // ---- Asymptomatic transitions ----
      add({iIa,    iIs},    {-1, +1});
      add({iIa,    iIso},   {-1, +1});
      add({iIa,    iRec},   {-1, +1});
      
      // ---- Symptomatic transitions ----
      add({iIs,    iIso},   {-1, +1});
      add({iIs,    iRec},   {-1, +1});
      
      // ---- Isolated transitions ----
      add({iIso,   iRec},   {-1, +1});        // 21 Iso   -> R    (isolated recovers)
      
      // ---- Waning immunity (three-stage) ----
      add({iR,     iR2},    {-1, +1});        // 22 R     -> R2
      add({iR2,    iR3},    {-1, +1});        // 23 R2    -> R3
      add({iR3,    iS},     {-1, +1});        // 24 R3    -> S
      
      // ---- Births ----
      add({iS},             {+1});            // 25 birth -> S
      
      // ---- Deaths (all human compartments) ----
      add({iS},     {-1});                    // 26
      add({iSprep}, {-1});                    // 27
      add({iSpost}, {-1});                    // 28
      add({iQs},    {-1});                    // 29
      add({iL},     {-1});                    // 30
      add({iQL},    {-1});                    // 31
      add({iIs},    {-1});                    // 32
      add({iIa},    {-1});                    // 33
      add({iIso},   {-1});                    // 34
      add({iR},     {-1});                    // 35
      add({iR2},    {-1});                    // 36
      add({iR3},    {-1});                    // 37 (cumul not mortal)
    }
  }
  
  // ---- Vector reactions (6, patch-global) ----
  auto add_g = [&](std::initializer_list<int> idxs,
                   std::initializer_list<int> dels) {
    Reaction r;
    r.idx.assign(idxs.begin(), idxs.end());
    r.delta.assign(dels.begin(), dels.end());
    rxn.push_back(std::move(r));
  };
  
  add_g({iSv, iEv, iCumulV}, {-1, +1, +1}); // V1 Sv -> Ev (vector infected by biting)
  add_g({iEv, iIv},           {-1, +1});    // V2 Ev -> Iv (extrinsic incubation complete)
  add_g({iSv},                {+1});        // V3 birth  -> Sv
  add_g({iSv},                {-1});        // V4 Sv death
  add_g({iEv},                {-1});        // V5 Ev death
  add_g({iIv},                {-1});        // V6 Iv death
  
  // Sanity check
  int expected = P * A * (RXN_PER_PA - (no_latent ? RXN_LATENT_DELTA : 0)) + RXN_VECTOR;
  if ((int)rxn.size() != expected)
    Rcpp::stop("build_reactions: expected %d reactions, got %d",
               expected, (int)rxn.size());
  
  return rxn;
}


// Compute force of infection lambda_{p,a} for every (p,a) pair.
// Combines:
//   (a) Human-to-human contact: beta_eff * seasonal * f_interv
//                               * sum_q Mp[p,q] * sum_b C[a,b] * Ieff_{q,b}/N_{q,b}
//   (b) Vector-to-human:        a_bite * b_h * (Iv/Nv) * beta_eff * seasonal * f_interv
//
// Infectious compartments: Is (full) and Ia (reduced by asymp_trans)
// Iso does NOT contribute (isolated = no contacts with humans or vectors)
static std::vector<double> compute_lambda_pa(
    const std::vector<double>& x,
    int P, int A,
    int offS,  int offSprep, int offSpost, int offQs,
    int offL,  int offQL,    int offIs,    int offIa,  int offIso,
    int offR,  int offR2,    int offR3,
    const NumericMatrix& C,
    const NumericMatrix& Mp,
    double beta_eff,
    double asymp_trans,
    double seasonal,
    double f_interv,
    double Iv, double Nv,
    double a_bite, double b_h
) {
  std::vector<double> lambda(P * A, 0.0);
  
  auto get = [&](int p, int a, int off) {
    return x[off + (p * A + a)];
  };
  
  // Vector-to-human component (same for all (p,a))
  double lambda_vec = 0.0;
  if (Nv > 0.0)
    lambda_vec = a_bite * b_h * (Iv / Nv) * beta_eff * seasonal * f_interv;
  
  for (int p = 0; p < P; ++p) {
    for (int a = 0; a < A; ++a) {
      
      double h2h = 0.0;
      
      for (int q = 0; q < P; ++q) {
        double mpq = Mp(p, q);
        if (mpq == 0.0) continue;
        
        for (int b = 0; b < A; ++b) {
          double S     = get(q, b, offS);
          double Sp    = get(q, b, offSprep);
          double Spost = get(q, b, offSpost);
          double Qs    = get(q, b, offQs);
          double L     = get(q, b, offL);
          double QL    = get(q, b, offQL);
          double Is    = get(q, b, offIs);
          double Ia    = get(q, b, offIa);
          double Iso   = get(q, b, offIso);
          double R     = get(q, b, offR);
          double R2    = get(q, b, offR2);
          double R3    = get(q, b, offR3);
          
          double N = S + Sp + Spost + Qs + L + QL + Is + Ia + Iso + R + R2 + R3;
          if (N <= 0.0) continue;
          
          // Iso does NOT transmit (isolated = no contacts)
          double Ieff = Is + asymp_trans * Ia;
          if (Ieff <= 0.0) continue;
          
          h2h += mpq * C(a, b) * (Ieff / N);
        }
      }
      
      lambda[p * A + a] = beta_eff * seasonal * f_interv * h2h + lambda_vec;
    }
  }
  
  return lambda;
}


// Compute force of infection on vectors FROM infected humans.
//   lambda_v = a_bite * b_v * (Ieff_total / Nh_total) * beta_eff * seasonal * f_interv
// Iso excluded (isolated humans have no vector contact).
static double compute_lambda_v(
    const std::vector<double>& x,
    int P, int A,
    int offS,  int offSprep, int offSpost, int offQs,
    int offL,  int offQL,    int offIs,    int offIa,  int offIso,
    int offR,  int offR2,    int offR3,
    double a_bite, double b_v, double asymp_trans,
    double beta_eff, double seasonal, double f_interv
) {
  auto get = [&](int p, int a, int off) {
    return x[off + (p * A + a)];
  };
  
  double Ieff_sum = 0.0;
  double Nh_total = 0.0;
  
  for (int p = 0; p < P; ++p) {
    for (int a = 0; a < A; ++a) {
      double S     = get(p, a, offS);
      double Sp    = get(p, a, offSprep);
      double Spost = get(p, a, offSpost);
      double Qs    = get(p, a, offQs);
      double L     = get(p, a, offL);
      double QL    = get(p, a, offQL);
      double Is    = get(p, a, offIs);
      double Ia    = get(p, a, offIa);
      double Iso   = get(p, a, offIso);
      double R     = get(p, a, offR);
      double R2    = get(p, a, offR2);
      double R3    = get(p, a, offR3);
      
      double N = S + Sp + Spost + Qs + L + QL + Is + Ia + Iso + R + R2 + R3;
      if (N <= 0.0) continue;
      
      Ieff_sum += Is + asymp_trans * Ia; // Iso excluded
      Nh_total += N;
    }
  }
  
  if (Nh_total <= 0.0) return 0.0;
  return a_bite * b_v * (Ieff_sum / Nh_total) * beta_eff * seasonal * f_interv;
}


// Compute propensities in SAME ORDER as build_reactions.
static void compute_propensities(
    std::vector<double>& a,
    double& totalRate,
    double& critRate,
    std::vector<char>& isCrit,
    const std::vector<double>& x,
    int P, int A,
    int offS,    int offSprep, int offSpost, int offQs,
    int offL,    int offQL,    int offIs,    int offIa,  int offIso,
    int offR,    int offR2,    int offR3,    int offCumul,
    int iSv, int iEv, int iIv,
    const std::vector<double>& lambda_pa,
    double lambda_v,
    double prep_start_rate_day,
    double prep_stop_rate_day,
    double quar_rate_day,
    double leave_quar_rate_day,
    double leave_quar_prep_rate_day,
    double pep_rate_day,
    double spost_waning_rate,
    double sigma,
    double gamma,
    double omega,
    double p_asym,
    double asymp_trans,
    double kappa,
    double iso_rate_day,
    double prep_effective_day,
    double birth_rate,
    double death_rate,
    double sigma_v,
    double mu_v,
    int no_latent,
    int Ncritical,
    const std::vector<Reaction>& rxn
) {
  int n_rxn = (int)rxn.size();
  a.assign(n_rxn, 0.0);
  isCrit.assign(n_rxn, 0);
  
  auto get = [&](int p, int aa, int off) {
    return x[off + (p * A + aa)];
  };
  
  totalRate = 0.0;
  critRate  = 0.0;
  int j     = 0;
  
  for (int p = 0; p < P; ++p) {
    for (int aa = 0; aa < A; ++aa) {
      
      double S     = get(p, aa, offS);
      double Sp    = get(p, aa, offSprep);
      double Spost = get(p, aa, offSpost);
      double Qs    = get(p, aa, offQs);
      double L     = get(p, aa, offL);
      double QL    = get(p, aa, offQL);
      double Is    = get(p, aa, offIs);
      double Ia    = get(p, aa, offIa);
      double Iso   = get(p, aa, offIso);
      double R     = get(p, aa, offR);
      double R2    = get(p, aa, offR2);
      double R3    = get(p, aa, offR3);
      double N     = S + Sp + Spost + Qs + L + QL + Is + Ia + Iso + R + R2 + R3;
      
      double lam = lambda_pa[p * A + aa];
      
      // ---- Infection ----
      // ---- Infection (mirrors build_reactions branch exactly) ----
      if (!no_latent) {
        a[j++] = lam * S;
        a[j++] = lam * (1.0 - prep_effective_day) * Sp;
      } else {
        // Total infection flux identical: lam*S split (1-p_asym):p_asym
        a[j++] = lam * (1.0 - p_asym) * S;
        a[j++] = lam * p_asym * S;
        a[j++] = lam * (1.0 - prep_effective_day) * (1.0 - p_asym) * Sp;
        a[j++] = lam * (1.0 - prep_effective_day) * p_asym * Sp;
      }
      
      // ---- PrEP transitions ----
      a[j++] = prep_start_rate_day * S;
      a[j++] = prep_stop_rate_day  * Sp;
      
      // ---- Quarantine (susceptible) ----
      a[j++] = quar_rate_day            * S;
      a[j++] = leave_quar_rate_day      * Qs;
      a[j++] = leave_quar_prep_rate_day * Qs;
      
      if (!no_latent) {
        a[j++] = quar_rate_day          * L;      // L -> QL
      }
      
      // ---- Sprep quarantine ----
      a[j++] = quar_rate_day            * Sp;
      
      if (!no_latent) {
        a[j++] = pep_rate_day      * L;
        a[j++] = pep_rate_day      * QL;
      }
      a[j++] = spost_waning_rate * Spost;
      
      if (!no_latent) {
        a[j++] = sigma * (1.0 - p_asym) * L;
        a[j++] = sigma * p_asym          * L;
        a[j++] = sigma * p_asym          * QL;
        a[j++] = sigma * (1.0 - p_asym) * QL;
      }
      
      // ---- Asymptomatic transitions ----
      a[j++] = kappa        * Ia;
      a[j++] = iso_rate_day * Ia;
      a[j++] = gamma        * Ia;
      
      // ---- Symptomatic transitions ----
      a[j++] = iso_rate_day * Is;
      a[j++] = gamma        * Is;
      
      // ---- Isolated transitions ----
      a[j++] = gamma        * Iso;                   // 21 Iso   -> R
      
      // ---- Waning immunity (three-stage) ----
      a[j++] = omega * R;                            // 22 R     -> R2
      a[j++] = omega * R2;                           // 23 R2    -> R3
      a[j++] = omega * R3;                           // 24 R3    -> S
      
      // ---- Births ----
      a[j++] = birth_rate * N;                       // 25 birth -> S
      
      // ---- Deaths (all human compartments) ----
      a[j++] = death_rate * S;                       // 26
      a[j++] = death_rate * Sp;                      // 27
      a[j++] = death_rate * Spost;                   // 28
      a[j++] = death_rate * Qs;                      // 29
      a[j++] = death_rate * L;                       // 30
      a[j++] = death_rate * QL;                      // 31
      a[j++] = death_rate * Is;                      // 32
      a[j++] = death_rate * Ia;                      // 33
      a[j++] = death_rate * Iso;                     // 34
      a[j++] = death_rate * R;                       // 35
      a[j++] = death_rate * R2;                      // 36
      a[j++] = death_rate * R3;                      // 37
    }
  }
  
  // ---- Vector propensities (6, matching vector reactions in build_reactions) ----
  double Sv = x[iSv];
  double Ev = x[iEv];
  double Iv = x[iIv];
  double Nv = Sv + Ev + Iv;
  
  a[j++] = lambda_v * Sv;  // V1 Sv -> Ev (vector infected)
  a[j++] = sigma_v  * Ev;  // V2 Ev -> Iv (extrinsic incubation)
  a[j++] = mu_v     * Nv;  // V3 birth -> Sv
  a[j++] = mu_v     * Sv;  // V4 Sv death
  a[j++] = mu_v     * Ev;  // V5 Ev death
  a[j++] = mu_v     * Iv;  // V6 Iv death
  
  // Clamp and accumulate totals
  for (int k = 0; k < n_rxn; k++) {
    double ak = a[k];
    if (!R_finite(ak) || ak < 0.0) ak = 0.0;
    a[k] = ak;
    totalRate += ak;
    if (ak > 0.0 && is_critical(x, rxn[k], Ncritical)) {
      isCrit[k] = 1;
      critRate  += ak;
    }
  }
}


// Cao tau-leaping step size selection
static double compute_cao_tau(
    const std::vector<double>& x,
    const std::vector<double>& a,
    const std::vector<char>& isCrit,
    const std::vector<Reaction>& rxn,
    double epsilon,
    const std::vector<double>& g,
    double maxtau
) {
  int M = (int)x.size();
  std::vector<double> mu(M, 0.0), sig(M, 0.0);
  
  for (int j = 0; j < (int)rxn.size(); ++j) {
    if (isCrit[j]) continue;
    double aj = a[j];
    if (aj <= 0.0) continue;
    const Reaction& r = rxn[j];
    for (size_t k = 0; k < r.idx.size(); ++k) {
      int i    = r.idx[k];
      int vij  = r.delta[k];
      mu[i]   += (double)vij * aj;
      sig[i]  += (double)(vij * vij) * aj;
    }
  }
  
  auto tau_bound = [&](double X, double gi, double mui, double sigi) {
    double xi = std::max(epsilon * X / gi, 1.0);
    double t1 = R_PosInf, t2 = R_PosInf;
    if (mui  != 0.0) t1 = xi / std::fabs(mui);
    if (sigi >  0.0) t2 = (xi * xi) / sigi;
    return std::min(t1, t2);
  };
  
  double tau = R_PosInf;
  for (int i = 0; i < M; ++i) {
    double gi = (g[i] > 0.0) ? g[i] : 1.0;
    double tb = tau_bound(x[i], gi, mu[i], sig[i]);
    if (tb < tau) tau = tb;
  }
  if (!R_finite(tau) || tau <= 0.0) tau = 1e-6;
  if (R_finite(maxtau)) tau = std::min(tau, maxtau);
  return tau;
}


}  // anonymous namespace

// [[Rcpp::export]]
List simulate_vectorborne_atl_cpp(List cfg) {
  
  // ---- dimensions ----
  int P       = as<int>(cfg["P"]);
  int A       = as<int>(cfg["A"]);
  int tf_days = as<int>(cfg["tf_days"]);

  // ---- Structural model flags ----
  int no_latent = cfg.containsElementNamed("no_latent") ? as<int>(cfg["no_latent"]) : 0;
  int sis_mode  = cfg.containsElementNamed("sis_mode")  ? as<int>(cfg["sis_mode"])  : 0;
  
  if (P <= 0) Rcpp::stop("P must be > 0");
  if (A <= 0) Rcpp::stop("A must be > 0");
  if (tf_days <= 0) Rcpp::stop("tf_days must be > 0");
  
  if (!cfg.containsElementNamed("C"))  Rcpp::stop("cfg missing 'C' (AxA contact matrix)");
  if (!cfg.containsElementNamed("Mp")) Rcpp::stop("cfg missing 'Mp' (PxP mixing matrix)");
  
  NumericMatrix C  = as<NumericMatrix>(cfg["C"]);
  NumericMatrix Mp = as<NumericMatrix>(cfg["Mp"]);
  
  if (C.nrow()  != A || C.ncol()  != A)
    Rcpp::stop("C must be %d x %d", A, A);
  if (Mp.nrow() != P || Mp.ncol() != P)
    Rcpp::stop("Mp must be %d x %d", P, P);
  
  // ---- daily driver vectors ----
  NumericVector beta_eff_daily        = as<NumericVector>(cfg["beta_eff_daily"]);
  NumericVector importation_daily     = as<NumericVector>(cfg["importation_daily"]);
  NumericVector seasonal_daily        = as<NumericVector>(cfg["seasonal_daily"]);
  NumericVector f_interv_daily        = as<NumericVector>(cfg["f_interv_daily"]);
  NumericVector iso_rate_daily        = as<NumericVector>(cfg["iso_rate_daily"]);
  NumericVector quar_rate_daily       = as<NumericVector>(cfg["quar_rate_daily"]);
  NumericVector leave_quar_daily      = as<NumericVector>(cfg["leave_quar_daily"]);
  NumericVector leave_quar_prep_daily = as<NumericVector>(cfg["leave_quar_prep_daily"]);
  NumericVector prep_start_daily      = as<NumericVector>(cfg["prep_start_daily"]);
  NumericVector prep_stop_daily       = as<NumericVector>(cfg["prep_stop_daily"]);
  NumericVector pep_rate_daily        = as<NumericVector>(cfg["pep_rate_daily"]);
  NumericVector prep_eff_daily        = as<NumericVector>(cfg["prep_eff_daily"]);
  // ---- Reactive intervention control (threshold-triggered) ----
  int    reactive_mode   = cfg.containsElementNamed("reactive_mode")   ? as<int>(cfg["reactive_mode"])      : 0;
  double on_threshold    = cfg.containsElementNamed("on_threshold")    ? as<double>(cfg["on_threshold"])    : 0.0;
  double off_threshold   = cfg.containsElementNamed("off_threshold")   ? as<double>(cfg["off_threshold"])   : 0.0;
  int    trigger_delay   = cfg.containsElementNamed("trigger_delay")   ? as<int>(cfg["trigger_delay"])      : 0;
  double reactive_contact_mult = cfg.containsElementNamed("reactive_contact_mult")
                                 ? as<double>(cfg["reactive_contact_mult"]) : 1.0;

  bool interv_active     = false;
  int  pending_countdown = -1;
  double prev_day_cases  = 0.0;
  NumericVector interv_active_daily(tf_days, 0.0);

  
  // ---- initial state ----
  // 13*P*A human entries + 4 vector entries
  NumericVector x0r   = as<NumericVector>(cfg["x0"]);
  int block   = P * A;
  int M_human = block * N_HUMAN_BLOCKS;
  int M       = M_human + 4;
  
  if ((int)x0r.size() != M)
    Rcpp::stop("x0 length mismatch: expected %d (= P*A*%d + 4), got %d",
               M, N_HUMAN_BLOCKS, (int)x0r.size());
  
  // Human offsets
  int offS     =  0 * block;
  int offSprep =  1 * block;
  int offSpost =  2 * block;
  int offQs    =  3 * block;
  int offL     =  4 * block;
  int offQL    =  5 * block;
  int offIs    =  6 * block;
  int offIa    =  7 * block;
  int offIso   =  8 * block;
  int offR     =  9 * block;
  int offR2    = 10 * block;
  int offR3    = 11 * block;
  int offCumul = 12 * block;
  
  // Vector offsets (appended after human blocks)
  int iSv     = M_human + 0;
  int iEv     = M_human + 1;
  int iIv     = M_human + 2;
  int iCumulV = M_human + 3;
  
  // ---- biological parameters ----
  List par = cfg["params"];
  double sigma             = as<double>(par["sigma"]);
  double gamma             = as<double>(par["gamma"]);
  double omega             = as<double>(par["omega"]);
  double p_asym            = as<double>(par["p_asym"]);
  double asymp_trans       = as<double>(par["asymp_trans"]);
  double kappa             = as<double>(par["kappa"]);
  double birth_rate        = as<double>(par["birth_rate"]);
  double death_rate        = as<double>(par["death_rate"]);
  double spost_waning_rate = as<double>(par["spost_waning_rate"]);
  
  // Vector parameters
  double sigma_v = as<double>(par["sigma_v"]);
  double mu_v    = as<double>(par["mu_v"]);
  double a_bite  = as<double>(par["a_bite"]);
  double b_h     = as<double>(par["b_h"]);
  double b_v     = as<double>(par["b_v"]);
  
  // ---- tau-leaping controls ----
  double epsilon        = cfg.containsElementNamed("epsilon")        ? as<double>(cfg["epsilon"])        : 0.03;
  int    Ncritical      = cfg.containsElementNamed("Ncritical")      ? as<int>(cfg["Ncritical"])          : 10;
  double exactThreshold = cfg.containsElementNamed("exactThreshold") ? as<double>(cfg["exactThreshold"])  : 10.0;
  double maxtau         = cfg.containsElementNamed("maxtau")         ? as<double>(cfg["maxtau"])          : R_PosInf;
  int    maxHalvings    = cfg.containsElementNamed("maxHalvings")    ? as<int>(cfg["maxHalvings"])        : 20;
  int    maxStepsPerDay = cfg.containsElementNamed("maxStepsPerDay") ? as<int>(cfg["maxStepsPerDay"])     : 10000000;
  
  std::vector<double> g(M, 1.0);
  if (cfg.containsElementNamed("changeBound")) {
    NumericVector gv = as<NumericVector>(cfg["changeBound"]);
    if ((int)gv.size() == M) for (int i = 0; i < M; ++i) g[i] = gv[i];
  }
  
  // ---- build reactions ----
  std::vector<Reaction> rxn = build_reactions(
    P, A, no_latent, sis_mode,
    offS, offSprep, offSpost, offQs,
    offL, offQL, offIs, offIa, offIso,
    offR, offR2, offR3, offCumul,
    iSv, iEv, iIv, iCumulV
  );
  int n_rxn = (int)rxn.size();
  
  // ---- state ----
  std::vector<double> x(M);
  for (int i = 0; i < M; ++i) x[i] = x0r[i];
  
  // ---- outputs ----
  NumericMatrix daily_cases_block(tf_days, block);
  NumericMatrix state_over_time(tf_days, M);
  NumericVector daily_cases_vector(tf_days, 0.0);
  NumericVector cumul_sum(tf_days + 1, 0.0);
  
  cumul_sum[0] = 0.0;
  for (int i = 0; i < block; ++i) cumul_sum[0] += x[offCumul + i];
  
  std::vector<double> prev_cumul(block);
  for (int i = 0; i < block; ++i) prev_cumul[i] = x[offCumul + i];
  double prev_cumul_v = x[iCumulV];
  
  // ---- scratch ----
  std::vector<double> a_prop;
  std::vector<char>   isCrit;
  double totalRate = 0.0, critRate = 0.0;
  
  auto drv = [&](const NumericVector& v, int day) -> double {
    return v[std::min(day, (int)v.size() - 1)];
  };
  
  // ---- main simulation loop ----
  for (int day = 0; day < tf_days; ++day) {
    
    double beta_eff            = drv(beta_eff_daily,        day);
    double import_rate         = drv(importation_daily,     day);  
    double seasonal            = drv(seasonal_daily,        day);
    double f_interv            = drv(f_interv_daily,        day);
    double iso_rate_day        = drv(iso_rate_daily,        day);
    double quar_rate_day       = drv(quar_rate_daily,       day);
    double leave_quar_day      = drv(leave_quar_daily,      day);
    double leave_quar_prep_day = drv(leave_quar_prep_daily, day);
    double prep_start_day      = drv(prep_start_daily,      day);
    double prep_stop_day       = drv(prep_stop_daily,       day);
    double pep_rate_day        = drv(pep_rate_daily,        day);
    double prep_eff_day        = drv(prep_eff_daily,        day);
    
    double t_in_day  = 0.0;
    int    step_count = 0;

    // ---- Reactive intervention state machine ----
    // Trigger signal is the PREVIOUS day's new cases, so the decision for
    // day d uses only information available at the end of day d-1.
    if (reactive_mode) {
      if (!interv_active) {
        if (pending_countdown < 0) {
          if (prev_day_cases >= on_threshold) pending_countdown = trigger_delay;
        } else if (pending_countdown > 0) {
          pending_countdown--;
        }
        if (pending_countdown == 0) { interv_active = true; pending_countdown = -1; }
      } else {
        if (prev_day_cases <= off_threshold) interv_active = false;
      }

      if (interv_active) {
        beta_eff *= reactive_contact_mult;
      } else {
        iso_rate_day  = 0.0;
        quar_rate_day = 0.0;
        pep_rate_day  = 0.0;
        prep_eff_day  = 0.0;
      }
    }
    interv_active_daily[day] = interv_active ? 1.0 : 0.0;

    
    if (import_rate > 0.0) {
      int n_imports = (int)R::rpois(import_rate);
      if (n_imports > 0) {
        double S_total = 0.0;
        for (int i = 0; i < block; ++i) S_total += x[offS + i];
        if (S_total > 0.0) {
          int remaining = n_imports;
          for (int i = 0; i < block; ++i) {
            int n_cell;
            if (i < block - 1) {
              double frac = x[offS + i] / S_total;
              n_cell = (int)R::rbinom(remaining, frac);
            } else {
              n_cell = remaining;
            }
            remaining -= n_cell;
            double actual = std::min((double)n_cell, x[offS + i]);
            x[offS + i] -= actual;
            if (no_latent) {
              // NO-LATENT: imported cases enter Is / Ia directly
              double n_asym = R::rbinom(actual, p_asym);
              x[offIa + i] += n_asym;
              x[offIs + i] += (actual - n_asym);
            } else {
              x[offL + i] += actual;
            }
            x[offCumul + i] += actual;
          }
        }
      }
    }
    
    while (t_in_day < 1.0) {
      
      if (++step_count > maxStepsPerDay) {
        Rcpp::warning("Day %d: step limit hit at t=%.8f", day, t_in_day);
        break;
      }
      
      double Iv_now = x[iIv];
      double Nv_now = x[iSv] + x[iEv] + x[iIv];
      
      // Combined FOI for humans (h2h + vector-to-human)
      std::vector<double> lambda_pa = compute_lambda_pa(
        x, P, A,
        offS, offSprep, offSpost, offQs,
        offL, offQL, offIs, offIa, offIso,
        offR, offR2, offR3,
        C, Mp, beta_eff, asymp_trans, seasonal, f_interv,
        Iv_now, Nv_now, a_bite, b_h
      );
      
      // FOI for vectors (from infectious humans)
      double lambda_v = compute_lambda_v(
        x, P, A,
        offS, offSprep, offSpost, offQs,
        offL, offQL, offIs, offIa, offIso,
        offR, offR2, offR3,
        a_bite, b_v, asymp_trans,
        beta_eff, seasonal, f_interv
      );
      
      compute_propensities(
        a_prop, totalRate, critRate, isCrit,
        x, P, A,
        offS, offSprep, offSpost, offQs,
        offL, offQL, offIs, offIa, offIso,
        offR, offR2, offR3, offCumul,
        iSv, iEv, iIv,
        lambda_pa, lambda_v,
        prep_start_day, prep_stop_day,
        quar_rate_day, leave_quar_day, leave_quar_prep_day,
        pep_rate_day, spost_waning_rate,
        sigma, gamma, omega,
        p_asym, asymp_trans, kappa,
        iso_rate_day, prep_eff_day,
        birth_rate, death_rate,
        sigma_v, mu_v,
        no_latent, Ncritical, rxn
      );
      
      if (totalRate <= 0.0) break;
      
      double remaining = 1.0 - t_in_day;
      double tau = compute_cao_tau(x, a_prop, isCrit, rxn, epsilon, g, maxtau);
      tau = std::min(tau, remaining);
      if (!(tau > 0.0)) tau = std::min(1e-6, remaining);
      
      // SSA fallback if tau too small
      if (tau < exactThreshold / totalRate) {
        double dt = rexp_rate(totalRate);
        if (dt > remaining) break;
        
        double u = R::runif(0.0, 1.0) * totalRate;
        double cdf = 0.0;
        int which = -1;
        for (int j2 = 0; j2 < n_rxn; ++j2) {
          cdf += a_prop[j2];
          if (u <= cdf) { which = j2; break; }
        }
        if (which >= 0) {
          const Reaction& r = rxn[which];
          bool ok = true;
          for (size_t k = 0; k < r.idx.size(); ++k)
            if (r.delta[k] < 0 && x[r.idx[k]] < 1.0) { ok = false; break; }
            if (ok) apply_reaction(x, r);
        }
        t_in_day += dt;
        continue;
      }
      
      // tau-leaping
      double tau2 = rexp_rate(critRate);
      double dt   = std::min({tau, tau2, remaining});
      
      bool accepted = false;
      std::vector<double> x0 = x;
      
      for (int halv = 0; halv <= maxHalvings; ++halv) {
        x = x0;
        bool neg = false;
        
        for (int j2 = 0; j2 < n_rxn; ++j2) {
          if (isCrit[j2]) continue;
          double mean = a_prop[j2] * dt;
          if (mean <= 0.0) continue;
          double k = R::rpois(mean);
          if (k <= 0.0) continue;
          
          const Reaction& r = rxn[j2];
          for (size_t kk = 0; kk < r.idx.size(); ++kk)
            if (r.delta[kk] < 0 && x[r.idx[kk]] < (-r.delta[kk]) * k)
            { neg = true; break; }
            if (neg) break;
            
            for (size_t kk = 0; kk < r.idx.size(); ++kk)
              x[r.idx[kk]] += (double)r.delta[kk] * k;
        }
        
        if (!neg)
          for (int i = 0; i < M; ++i)
            if (x[i] < 0.0) { neg = true; break; }
            
            if (!neg) { accepted = true; break; }
            dt *= 0.5;
            if (!(dt > 0.0)) break;
      }
      
      // fallback to SSA if halving exhausted
      if (!accepted) {
        double dt2 = rexp_rate(totalRate);
        if (dt2 > remaining) break;
        
        double u = R::runif(0.0, 1.0) * totalRate;
        double cdf = 0.0;
        int which = -1;
        for (int j2 = 0; j2 < n_rxn; ++j2) {
          cdf += a_prop[j2];
          if (u <= cdf) { which = j2; break; }
        }
        if (which >= 0) {
          const Reaction& r = rxn[which];
          bool ok = true;
          for (size_t k = 0; k < r.idx.size(); ++k)
            if (r.delta[k] < 0 && x[r.idx[k]] < 1.0) { ok = false; break; }
            if (ok) apply_reaction(x, r);
        }
        t_in_day += dt2;
        continue;
      }
      
      t_in_day += dt;
      
      // execute ONE critical event if it fires before tau
      if (critRate > 0.0 && dt == tau2 && t_in_day < 1.0) {
        double u = R::runif(0.0, 1.0) * critRate;
        double cdf = 0.0;
        int which = -1;
        for (int j2 = 0; j2 < n_rxn; ++j2) {
          if (!isCrit[j2]) continue;
          cdf += a_prop[j2];
          if (u <= cdf) { which = j2; break; }
        }
        if (which >= 0) {
          const Reaction& r = rxn[which];
          bool ok = true;
          for (size_t k = 0; k < r.idx.size(); ++k)
            if (r.delta[k] < 0 && x[r.idx[k]] < 1.0) { ok = false; break; }
            if (ok) apply_reaction(x, r);
        }
      }
      
    } // end within-day loop
    
    // record daily new cases
    double cum_after = 0.0;
    for (int i = 0; i < block; ++i) {
      double new_cases = x[offCumul + i] - prev_cumul[i];
      if (!R_finite(new_cases) || new_cases < 0.0) new_cases = 0.0;
      daily_cases_block(day, i) = new_cases;
      prev_cumul[i] = x[offCumul + i];
      cum_after    += x[offCumul + i];
    }
    cumul_sum[day + 1] = cum_after;

    prev_day_cases = 0.0;
    for (int i = 0; i < block; ++i) prev_day_cases += daily_cases_block(day, i);

    
    // record daily new vector infections
    double new_vec = x[iCumulV] - prev_cumul_v;
    if (!R_finite(new_vec) || new_vec < 0.0) new_vec = 0.0;
    daily_cases_vector[day] = new_vec;
    prev_cumul_v = x[iCumulV];
    
    // save full state for this day
    for (int i = 0; i < M; ++i) state_over_time(day, i) = x[i];
    
  } // end day loop
  
  NumericVector x_final(M);
  for (int i = 0; i < M; ++i) x_final[i] = x[i];
  x_final.attr("names") = x0r.attr("names");
  
  return List::create(
    _["daily_new_cases_block"]  = daily_cases_block,
    _["daily_new_cases_vector"] = daily_cases_vector,
    _["cumul_infections"]       = cumul_sum,
    _["state_over_time"]        = state_over_time,
    _["intervention_active"]    = interv_active_daily,
    _["final_state"]            = x_final
  );
}