// mixbayes_gibbs.cpp
//
// Gibbs sampler with data augmentation for the two-component Bayesian mixture
// regression models behind glmMixBayes() and survregMixBayes()
// (Gutman et al., 2016).
//
// Model
//   z_n in {1, 2}            latent match indicator (1 = correct, 2 = mismatch)
//   P(z_n = 1) = theta                       Path A: theta ~ Beta(a, b)
//              = inv_logit(Z_n' gamma)       Path B: gamma ~ Normal
//   y_n | z_n = j ~ f_j(y_n | X_n' beta_j, dispersion_j)
//   "safe" records are fixed at z_n = 1 and do not enter the theta/gamma
//   likelihood (they still inform beta_1).
//
// One sweep of the sampler (component parameters given z, then the mixing
// weight, then z given the new parameter values)
//   beta_j, dispersion_j  given z (see below)
//   theta               collapsed update: slice sampling of logit(theta)
//                       with z integrated out, followed by an exact draw of
//                       z | theta (a blocked update of (theta, z))  [Path A]
//   gamma               Laplace independence MH given z    [Path B]
//   z_n                 exact Bernoulli draw
//   beta_j              exact Normal draw                  [gaussian]
//                       Laplace independence MH            [other families]
//   sigma / phi / shape slice sampling on the log scale
//   weibull scale       log(scale_j) is appended to the beta_j block so the
//                       ridge between scale and the intercept is handled by
//                       the joint proposal
//   all continuous      joint random-walk and independence Metropolis moves
//   parameters          with z integrated out (partially collapsed Gibbs),
//                       proposals estimated during the burn-in and fixed
//                       while draws are stored (see joint_move())
//
// Before the main chain, a few short pilot chains are run from different
// starting configurations (default initial values, mismatch component at the
// marginal outcome distribution, random allocations) and the main chain
// continues from the pilot with the highest average log posterior. This guards
// against the Gibbs sampler settling in a minor mode of the mixture posterior.
// When the R code asks for it (majority convention; see 'Label switching' in
// ?glmMixBayes), the labels of that starting state are exchanged if component
// 2 holds the majority of the records.
//
// The Laplace proposal is a multivariate t centred at the mode of the block's
// conditional posterior (found by damped Newton) with scale matrix equal to
// the inverse negative Hessian. With normal priors on the coefficients and a
// gamma / exponential / lognormal prior on the Weibull scale every block target
// is strictly log-concave, so the mode is unique; the Newton search always
// starts from the mode found at the previous update (never from the current
// draw), so the proposal depends on the conditioning variables only and the
// kernel is a valid independence Metropolis-Hastings step.

#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <chrono>
#include <climits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::vec;
using arma::mat;
using arma::uvec;
using arma::ivec;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();
const double LOG_2PI = 1.837877066409345483560659472811;
const double PROPOSAL_DF = 8.0;   // degrees of freedom of the t proposal

enum Family { GAUSSIAN = 1, POISSON = 2, BINOMIAL = 3, GAMMA = 4,
              SURV_GAMMA = 5, SURV_WEIBULL = 6 };

inline double log1p_exp(double x) {
  if (x > 35.0) return x + std::exp(-x);
  if (x < -35.0) return std::exp(x);
  return std::log1p(std::exp(x));
}

inline double inv_logit(double x) {
  if (x >= 0.0) { double e = std::exp(-x); return 1.0 / (1.0 + e); }
  double e = std::exp(x); return e / (1.0 + e);
}

// log S(t) = log Q(s, x) of a right-censored Gamma(shape s, mean exp(eta))
// record, with x = s * t * exp(-eta) computed on the log scale: x = Inf when
// exp(eta) underflows (log S = -Inf) and x = 0 when it overflows (log S = 0),
// so the result is never NaN.
inline double gamma_log_surv(double s, double logt, double eta) {
  return R::pgamma(std::exp(std::log(s) + logt - eta), s, 1.0, 0, 1);
}

// Derivatives of log Q(s, x) with respect to eta for x = s * t * exp(-eta):
// g = r and h = -r (s - x + r), where r = D / Q is the ratio of the
// standardised gamma density term D = x^s exp(-x) / Gamma(s) to Q. For large
// x, r and s - x + r are computed from the asymptotic expansion
// Q(s, x) ~ x^(s-1) exp(-x) / Gamma(s) * sum_k (s-1)...(s-k) / x^k, which
// avoids the cancellation of two numbers of size x.
inline void gamma_log_surv_derivs(double s, double x, double logQ, double& g, double& h) {
  if (!std::isfinite(x)) { g = 0.0; h = 0.0; return; }
  double r, c;   // c = s - x + r
  if (x > std::max(100.0, 10.0 * s)) {
    // xtail = x * sum_{k >= 1} (s-1)...(s-k) / x^k, so that r = x^2 / (x + xtail)
    double term = s - 1.0, xtail = term;
    for (int k = 2; k < 40; ++k) {
      term *= (s - k) / x;
      xtail += term;
      if (std::fabs(term) <= 1e-17 * std::fabs(xtail) || term == 0.0) break;
    }
    r = x / (1.0 + xtail / x);
    c = s - xtail / (1.0 + xtail / x);
  } else {
    r = std::exp(s * std::log(x) - x - std::lgamma(s) - logQ);
    c = s - x + r;
  }
  g = r;
  h = -r * c;
  if (!(h < 0.0)) h = 0.0;   // log S is concave in eta; guard rounding
}

// ---------------------------------------------------------------------------
// Priors on positive scalar parameters (sigma, phi, shape, scale)
// ---------------------------------------------------------------------------
struct ScalarPrior {
  // 1 normal(mu, sd)  2 cauchy(loc, scale)  3 gamma(shape, rate)  4 exponential(rate)
  // 5 lognormal(meanlog, sdlog)  6 student_t(nu, mu, sd)  7 inv_gamma(shape, scale)
  // Distributions with support on the whole real line (1, 2, 6) are truncated at 0.
  int dist;
  double p1, p2, p3;

  // unnormalised log density on x > 0
  double logd(double x) const {
    if (!(x > 0.0) || !std::isfinite(x)) return NEG_INF;
    switch (dist) {
      case 1: { double u = (x - p1) / p2; return -0.5 * u * u; }
      case 2: { double u = (x - p1) / p2; return -std::log1p(u * u); }
      case 3: return (p1 - 1.0) * std::log(x) - p2 * x;
      case 4: return -p1 * x;
      case 5: { double u = (std::log(x) - p1) / p2; return -0.5 * u * u - std::log(x); }
      case 6: { double u = (x - p2) / p3; return -0.5 * (p1 + 1.0) * std::log1p(u * u / p1); }
      case 7: return -(p1 + 1.0) * std::log(x) - p2 / x;
      default: return NEG_INF;
    }
  }
  // log density of u = log(x), Jacobian included
  double logd_log(double u) const { return logd(std::exp(u)) + u; }
  // first and second derivatives of logd_log(u); closed form for the
  // log-concave families used for the Weibull scale, central differences otherwise
  void dlogd_log(double u, double& d1, double& d2) const {
    if (dist == 3) { double x = std::exp(u); d1 = p1 - p2 * x; d2 = -p2 * x; return; }
    if (dist == 5) { d1 = -(u - p1) / (p2 * p2); d2 = -1.0 / (p2 * p2); return; }
    const double eps = 1e-5;
    double f0 = logd_log(u), fp = logd_log(u + eps), fm = logd_log(u - eps);
    d1 = (fp - fm) / (2.0 * eps);
    d2 = (fp - 2.0 * f0 + fm) / (eps * eps);
  }
};

ScalarPrior scalar_prior(const List& priors, const std::string& name) {
  NumericVector v = priors[name];
  if (v.size() < 3) stop("Scalar prior '%s' must be c(code, p1, p2[, p3]).", name.c_str());
  ScalarPrior p;
  p.dist = static_cast<int>(v[0]);
  p.p1 = v[1];
  p.p2 = v[2];
  p.p3 = (v.size() > 3) ? v[3] : 0.0;
  if (p.dist < 1 || p.dist > 7) stop("Unknown scalar prior code for '%s'.", name.c_str());
  return p;
}

// ---------------------------------------------------------------------------
// Likelihood terms
// ---------------------------------------------------------------------------
struct Model {
  int family;
  mat X;          // N x K design matrix
  vec y;          // outcome (survival time for survival families)
  vec logy;       // log(y) (gamma / weibull)
  vec lfact;      // lgamma(y + 1) (poisson)
  ivec event;     // 1 = event, 0 = right-censored (survival only)
  int N, K;

  bool weibull() const { return family == SURV_WEIBULL; }

  // log-likelihood of observation n given its full linear predictor eta
  // (for weibull, eta already includes log(scale)) and scalar parameter s
  // (sigma / phi / shape).
  double loglik_one(int n, double eta, double s) const {
    switch (family) {
      case GAUSSIAN: {
        double r = (y[n] - eta) / s;
        return -std::log(s) - 0.5 * LOG_2PI - 0.5 * r * r;
      }
      case POISSON:
        return y[n] * eta - std::exp(eta) - lfact[n];
      case BINOMIAL:
        return y[n] * eta - log1p_exp(eta);
      case GAMMA:
        return gamma_event(n, eta, s);
      case SURV_GAMMA:
        if (event[n] == 1) return gamma_event(n, eta, s);
        return gamma_log_surv(s, logy[n], eta);   // log S(t)
      case SURV_WEIBULL: {
        double w = std::exp(s * (logy[n] - eta));
        double ll = -w;
        if (event[n] == 1) ll += std::log(s) + (s - 1.0) * logy[n] - s * eta;
        return ll;
      }
    }
    return NEG_INF;
  }

  // log-likelihood plus first and second derivatives w.r.t. eta
  void terms_one(int n, double eta, double s, double& ll, double& g, double& h) const {
    switch (family) {
      case GAUSSIAN: {
        double r = (y[n] - eta) / s;
        ll = -std::log(s) - 0.5 * LOG_2PI - 0.5 * r * r;
        g = r / s;
        h = -1.0 / (s * s);
        return;
      }
      case POISSON: {
        double mu = std::exp(eta);
        ll = y[n] * eta - mu - lfact[n];
        g = y[n] - mu;
        h = -mu;
        return;
      }
      case BINOMIAL: {
        // one exponential serves both the probability and log(1 + e^eta)
        double p, lse;
        if (eta >= 0.0) { double e = std::exp(-eta); p = 1.0 / (1.0 + e); lse = eta + std::log1p(e); }
        else            { double e = std::exp(eta);  p = e / (1.0 + e);   lse = std::log1p(e); }
        ll = y[n] * eta - lse;
        g = y[n] - p;
        h = -p * (1.0 - p);
        return;
      }
      case GAMMA:
      case SURV_GAMMA: {
        if (family == GAMMA || event[n] == 1) {
          double e = y[n] * std::exp(-eta);
          ll = gamma_event(n, eta, s);
          g = s * (e - 1.0);
          h = -s * e;
          return;
        }
        // right-censored: log S = log Q(s, x) with x = s * t * exp(-eta)
        double x = std::exp(std::log(s) + logy[n] - eta);
        ll = R::pgamma(x, s, 1.0, 0, 1);
        gamma_log_surv_derivs(s, x, ll, g, h);
        return;
      }
      case SURV_WEIBULL: {
        double w = std::exp(s * (logy[n] - eta));
        ll = -w;
        g = s * w;
        h = -s * s * w;
        if (event[n] == 1) {
          ll += std::log(s) + (s - 1.0) * logy[n] - s * eta;
          g -= s;
        }
        return;
      }
    }
    ll = NEG_INF; g = 0.0; h = 0.0;
  }

  // Gamma(shape = s, rate = s * exp(-eta)) log density (mean = exp(eta))
  double gamma_event(int n, double eta, double s) const {
    return s * std::log(s) - std::lgamma(s) + (s - 1.0) * logy[n]
           - s * eta - s * y[n] * std::exp(-eta);
  }
};

// ---------------------------------------------------------------------------
// Block prior: independent normals on the first `ng` coordinates and an
// optional scalar prior on coordinate `ng` (the weibull log-scale).
// ---------------------------------------------------------------------------
struct BlockPrior {
  vec mean, sd;
  int ng;
  const ScalarPrior* sp;
  // prior mean of the full block (log-scale coordinate, if any, starts at 0)
  vec mean_full(int d) const { vec m(d, arma::fill::zeros); m.head(ng) = mean; return m; }
};

// Log conditional posterior of a coefficient block (and derivatives) given
// the observations `idx`, design `A` (rows aligned with idx) and scalar s.
struct BlockTarget {
  const Model& m;
  const mat& A;
  const uvec& idx;
  double s;
  const BlockPrior& pr;
  mutable vec g, h;

  double operator()(const vec& b, vec* grad, mat* hess) const {
    const int n = idx.n_elem;
    vec eta = A * b;
    double lp = 0.0;
    if (grad) { g.set_size(n); h.set_size(n); }
    for (int i = 0; i < n; ++i) {
      double ll;
      if (grad) {
        double gi, hi;
        m.terms_one(idx[i], eta[i], s, ll, gi, hi);
        g[i] = gi; h[i] = hi;
      } else {
        ll = m.loglik_one(idx[i], eta[i], s);
      }
      lp += ll;
    }
    for (int i = 0; i < pr.ng; ++i) {
      double u = (b[i] - pr.mean[i]) / pr.sd[i];
      lp -= 0.5 * u * u;
    }
    if (pr.sp) lp += pr.sp->logd_log(b[pr.ng]);
    if (!std::isfinite(lp)) return NEG_INF;

    if (grad) {
      *grad = A.t() * g;
      *hess = A.t() * (A.each_col() % h);
      for (int i = 0; i < pr.ng; ++i) {
        double v = pr.sd[i] * pr.sd[i];
        (*grad)[i] -= (b[i] - pr.mean[i]) / v;
        (*hess)(i, i) -= 1.0 / v;
      }
      if (pr.sp) {   // scalar prior on the log-scale coordinate
        double d1, d2;
        pr.sp->dlogd_log(b[pr.ng], d1, d2);
        (*grad)[pr.ng] += d1;
        (*hess)(pr.ng, pr.ng) += d2;
      }
    }
    return lp;
  }
};

// Triangular solves after a successful Cholesky factorisation: no condition
// number estimate and no approximate (SVD) fallback.
inline vec solve_upper(const mat& R, const vec& b) {
  return arma::solve(arma::trimatu(R), b, arma::solve_opts::fast);
}
inline vec solve_lower(const mat& L, const vec& b) {
  return arma::solve(arma::trimatl(L), b, arma::solve_opts::fast);
}

// Upper Cholesky factor of a (nearly) positive definite matrix, P = R' R.
// The matrix is symmetrised and equilibrated (unit diagonal) first, so that a
// ridge added when the factorisation fails is relative to the scale of each
// direction rather than to the largest diagonal element.
bool chol_pd(mat P, mat& R) {
  P = 0.5 * (P + P.t());
  const vec d = P.diag();
  const bool equil = d.is_finite() && arma::all(d > 0.0);
  const vec s = equil ? vec(1.0 / arma::sqrt(d)) : vec(P.n_rows, arma::fill::ones);
  mat Q = equil ? mat(P % (s * s.t())) : P;
  const double base = equil ? 1e-8 : 1e-8 * std::max(1e-12, arma::abs(d).max());
  double ridge = 0.0;
  for (int attempt = 0; attempt < 12; ++attempt) {
    if (arma::chol(R, Q)) {
      if (equil) R.each_row() /= s.t();   // P = D^-1 Q D^-1 with D = diag(s)
      return true;
    }
    ridge = (ridge == 0.0) ? base : ridge * 10.0;
    Q.diag() += ridge;
  }
  return false;
}

// Damped Newton search for the mode of a strictly log-concave block target.
// Convergence is declared when the squared Newton decrement g' (-H)^-1 g,
// which approximates twice the log-density gap to the mode and does not
// depend on the scale of the coefficients, is negligible. On success `b` holds
// the mode and `P` the (ridged) negative Hessian; returns false when the
// search does not converge, so that the caller can reject.
bool find_mode(const BlockTarget& f, vec& b, mat& P, double& lp) {
  vec grad, gradn; mat H, Hn;
  lp = f(b, &grad, &H);
  if (!std::isfinite(lp)) return false;
  bool converged = false;
  for (int it = 0; it < 100; ++it) {
    mat R;
    if (!chol_pd(-H, R)) return false;
    vec step = solve_upper(R, solve_lower(R.t(), grad));
    const double decrement = arma::dot(grad, step);
    if (!std::isfinite(decrement)) return false;
    if (decrement < 1e-10) { converged = true; break; }
    // Take the full Newton step, halving it while the objective does not
    // improve. Derivatives are evaluated together with the objective so that
    // an accepted step costs a single evaluation.
    double t = 1.0;
    bool accepted = false;
    vec bn;
    for (int bt = 0; bt < 40; ++bt) {
      bn = b + t * step;
      double lpn = f(bn, &gradn, &Hn);
      if (std::isfinite(lpn) && lpn >= lp - 1e-12 * std::fabs(lp)) {
        b = bn; lp = lpn; grad = gradn; H = Hn; accepted = true;
        break;
      }
      t *= 0.5;
    }
    if (!accepted) break;
  }
  P = -H;
  return converged && std::isfinite(lp);
}

// Laplace independence Metropolis-Hastings update of one block.
// `start` is the Newton starting point: the mode found at the previous update
// (the initial value at the first update). The block target is strictly
// log-concave, so the mode that is found does not depend on the starting
// point; the proposal is then a function of the conditioning variables only
// and the kernel is a valid independence sampler. A failed mode search counts
// as a rejection.
bool laplace_mh(const BlockTarget& f, vec& b, vec& start) {
  const int d = b.n_elem;
  const double nu = PROPOSAL_DF;
  vec mode = start;
  if (!mode.is_finite() || static_cast<int>(mode.n_elem) != d) mode = f.pr.mean_full(d);
  mat P; double lp_mode;
  if (!find_mode(f, mode, P, lp_mode)) {
    mode = f.pr.mean_full(d);
    if (!find_mode(f, mode, P, lp_mode)) return false;
  }
  start = mode;
  mat R;
  if (!chol_pd(P, R)) return false;

  vec eps(d);
  for (int i = 0; i < d; ++i) eps[i] = norm_rand();
  double w = R::rchisq(nu);
  vec prop = mode + solve_upper(R, eps) * std::sqrt(nu / w);

  auto logq = [&](const vec& x) {
    vec Rd = R * (x - mode);
    return -0.5 * (nu + d) * std::log1p(arma::dot(Rd, Rd) / nu);
  };

  double lp_prop = f(prop, nullptr, nullptr);
  if (!std::isfinite(lp_prop)) return false;
  double lp_cur = f(b, nullptr, nullptr);
  double la = (lp_prop - logq(prop)) - (lp_cur - logq(b));
  if (!std::isfinite(lp_cur) || std::log(unif_rand()) < la) { b = prop; return true; }
  return false;
}

// Exact conjugate draw of beta for the gaussian family.
void gaussian_beta_draw(const Model& m, const uvec& idx, double sigma,
                        const BlockPrior& pr, vec& b) {
  const int K = m.K;
  mat P(K, K, arma::fill::zeros);
  vec r(K, arma::fill::zeros);
  if (idx.n_elem > 0) {
    mat A = m.X.rows(idx);
    vec yy = m.y.elem(idx);
    P = A.t() * A / (sigma * sigma);
    r = A.t() * yy / (sigma * sigma);
  }
  for (int i = 0; i < K; ++i) {
    double v = pr.sd[i] * pr.sd[i];
    P(i, i) += 1.0 / v;
    r[i] += pr.mean[i] / v;
  }
  mat R;
  if (!chol_pd(P, R)) stop("Cholesky factorisation failed in the gaussian beta update.");
  vec mu = solve_upper(R, solve_lower(R.t(), r));
  vec eps(K);
  for (int i = 0; i < K; ++i) eps[i] = norm_rand();
  b = mu + solve_upper(R, eps);
}

// Univariate slice sampler with stepping out and shrinkage (Neal, 2003).
template <class F>
double slice_sample(double x0, F&& logf, double w = 1.0, int max_steps = 50) {
  double f0 = logf(x0);
  if (!std::isfinite(f0)) return x0;
  double level = f0 - R::exp_rand();
  double L = x0 - w * unif_rand(), U = L + w;
  int j = static_cast<int>(std::floor(max_steps * unif_rand()));
  int k = max_steps - 1 - j;
  while (j > 0 && logf(L) > level) { L -= w; --j; }
  while (k > 0 && logf(U) > level) { U += w; --k; }
  for (int it = 0; it < 500; ++it) {
    double x1 = L + (U - L) * unif_rand();
    if (logf(x1) > level) return x1;
    if (x1 < x0) L = x1; else U = x1;
  }
  return x0;
}

// Sum of log-likelihood terms over `idx` for a given scalar parameter.
double component_loglik(const Model& m, const uvec& idx, const vec& eta, double s) {
  double lp = 0.0;
  for (arma::uword i = 0; i < idx.n_elem; ++i) lp += m.loglik_one(idx[i], eta[i], s);
  return lp;
}

// Sample SD of v, or the reference scale `ref` of the outcome when v is too
// short or (numerically) constant. Relative to `ref`, so that the starting
// values do not depend on the units of the outcome.
double sample_sd(const vec& v, double ref) {
  if (v.n_elem < 2) return ref;
  double s = arma::stddev(v);
  return (std::isfinite(s) && s > 1e-8 * ref) ? s : ref;
}

// Running log of the mean of exp(x) over the values added (log-sum-exp form);
// NaN once a NaN was added, NA_REAL when nothing was added.
struct LogMeanExp {
  double m = NEG_INF, s = 0.0;
  long n = 0;
  bool nan = false;
  void add(double x) {
    ++n;
    if (std::isnan(x)) { nan = true; return; }
    if (x == NEG_INF) return;
    if (x > m) { s = s * std::exp(m - x) + 1.0; m = x; } else s += std::exp(x - m);
  }
  double value() const {
    if (n == 0) return NA_REAL;
    if (nan) return R_NaN;
    return (m == NEG_INF) ? NEG_INF : m + std::log(s / static_cast<double>(n));
  }
};

// Reference scale of the outcome: its SD, else its mean absolute value, else 1.
double outcome_scale(const vec& y) {
  double s = (y.n_elem > 1) ? arma::stddev(y) : 0.0;
  if (std::isfinite(s) && s > 0.0) return s;
  double a = y.n_elem > 0 ? arma::mean(arma::abs(y)) : 0.0;
  return (std::isfinite(a) && a > 0.0) ? a : 1.0;
}

} // namespace

// ---------------------------------------------------------------------------
// Sampler entry point
// ---------------------------------------------------------------------------
// `pre_orient`: how the starting state of the main chain is oriented (see the
// orientation block below): 0 not at all (safe matches identify the labels);
// 1 so that component 1 holds at least half of the records (the majority
// convention, used by the R code when neither safe matches nor the prior on
// the match probability identify the labels; TRUE from R); 2 towards the
// labelling with the larger posterior probability (used when the prior on the
// match probability identifies the labels). Ignored with user starting values
// and with safe matches. `orient_tol` (majority convention only): the labels
// are exchanged when the estimated log ratio of the posterior probabilities
// of the exchanged and the original labelling is at least -orient_tol; below,
// check chains decide (Inf always exchanges them).
//' @noRd
// [[Rcpp::export]]
List mixbayes_gibbs_cpp(std::string family, const arma::mat& X, const arma::vec& y,
                        const arma::ivec& event, const arma::mat& Z, const arma::ivec& safe,
                        List priors, List init,
                        int n_iter, int n_burnin, int thin,
                        int n_pilot, int pilot_iter, bool collapse_theta, bool verbose,
                        int pre_orient = 0, double orient_tol = 2.0) {

  // ---- settings and the R output matrix of allocations ----------------------
  // The largest R allocation is made first, before any C++ object holds
  // memory, so that a failed allocation (an R error) cannot leak it.
  if (n_iter < 1 || n_burnin < 0 || n_burnin >= n_iter || thin < 1 || n_iter == INT_MAX)
    stop("Invalid iteration settings: need 0 <= burnin < iterations and thin >= 1.");
  if (n_pilot < 0 || pilot_iter < 1 || pilot_iter == INT_MAX) stop("Invalid pilot-chain settings.");
  if (std::isnan(orient_tol) || orient_tol < 0.0) stop("`orient_tol` must be non-negative.");
  if (pre_orient < 0 || pre_orient > 2) stop("`pre_orient` must be 0, 1 or 2.");
  const int S = (n_iter - n_burnin) / thin;
  if (S < 1) stop("No draws would be stored: increase `iterations` or reduce `thin`.");
  if (X.n_rows < 1) stop("`X` must have at least one row.");
  IntegerMatrix z_s(S, static_cast<int>(X.n_rows));

  // ---- model -------------------------------------------------------------
  Model m;
  if (family == "gaussian")            m.family = GAUSSIAN;
  else if (family == "poisson")        m.family = POISSON;
  else if (family == "binomial")       m.family = BINOMIAL;
  else if (family == "gamma")          m.family = GAMMA;
  else if (family == "surv_gamma")     m.family = SURV_GAMMA;
  else if (family == "surv_weibull")   m.family = SURV_WEIBULL;
  else stop("Unknown family '%s'.", family.c_str());

  const int N = X.n_rows, K = X.n_cols;
  if (K < 1) stop("`X` must have at least one column.");
  if (static_cast<int>(y.n_elem) != N) stop("`y` must have length nrow(X).");
  if (static_cast<int>(safe.n_elem) != N) stop("`safe.matches` must have length nrow(X).");
  if (!X.is_finite() || !y.is_finite()) stop("`X` and `y` must not contain missing or infinite values.");
  if (Z.n_cols > 0 && !Z.is_finite()) stop("`Z` must not contain missing or infinite values.");
  m.X = X; m.y = y; m.N = N; m.K = K;
  const bool survival = (m.family == SURV_GAMMA || m.family == SURV_WEIBULL);
  if (survival) {
    if (static_cast<int>(event.n_elem) != N) stop("`event` must have length nrow(X).");
    m.event = event;
  }
  if (m.family == GAMMA || survival) {
    if (y.min() <= 0.0) stop("Outcomes must be strictly positive for this family.");
    m.logy = arma::log(y);
  }
  if (m.family == POISSON) {
    m.lfact.set_size(N);
    for (int n = 0; n < N; ++n) m.lfact[n] = std::lgamma(y[n] + 1.0);
  }
  const bool weibull = m.weibull();
  const bool use_logistic = Z.n_cols > 0;
  const int M = use_logistic ? Z.n_cols : 0;
  if (use_logistic && static_cast<int>(Z.n_rows) != N) stop("`Z` must have nrow(X) rows.");
  // whether column 1 of X is an intercept (set by the R code; used only for
  // starting values of the mismatch component)
  const bool has_intercept = priors.containsElementNamed("intercept") ? as<bool>(priors["intercept"]) : true;
  const double y_scale = outcome_scale(y);

  const uvec nonsafe = arma::find(safe == 0);
  const uvec all_idx = arma::regspace<uvec>(0, N - 1);

  // ---- priors ------------------------------------------------------------
  BlockPrior pb1, pb2, pg;
  pb1.mean = as<vec>(priors["beta1_mean"]); pb1.sd = as<vec>(priors["beta1_sd"]);
  pb2.mean = as<vec>(priors["beta2_mean"]); pb2.sd = as<vec>(priors["beta2_sd"]);
  pb1.ng = pb2.ng = K; pb1.sp = pb2.sp = nullptr;
  if (static_cast<int>(pb1.mean.n_elem) != K || static_cast<int>(pb1.sd.n_elem) != K ||
      static_cast<int>(pb2.mean.n_elem) != K || static_cast<int>(pb2.sd.n_elem) != K)
    stop("Coefficient prior vectors must have length ncol(X).");

  double theta_a = 1.0, theta_b = 1.0;
  if (use_logistic) {
    pg.mean = as<vec>(priors["gamma_mean"]); pg.sd = as<vec>(priors["gamma_sd"]);
    pg.ng = M; pg.sp = nullptr;
    if (static_cast<int>(pg.mean.n_elem) != M || static_cast<int>(pg.sd.n_elem) != M)
      stop("gamma prior vectors must have length ncol(Z).");
  } else {
    NumericVector th = priors["theta"];
    theta_a = th[0]; theta_b = th[1];
  }

  // scalar priors: disp = sigma (gaussian), phi (gamma), shape (weibull)
  ScalarPrior sp1, sp2, sc1, sc2;
  const bool has_disp = (m.family == GAUSSIAN || m.family == GAMMA || survival);
  if (has_disp) { sp1 = scalar_prior(priors, "disp1"); sp2 = scalar_prior(priors, "disp2"); }
  if (weibull) {
    sc1 = scalar_prior(priors, "scale1"); sc2 = scalar_prior(priors, "scale2");
    pb1.sp = &sc1; pb2.sp = &sc2;
  }

  // ---- state ---------------------------------------------------------------
  const int D = K + (weibull ? 1 : 0);      // size of the beta block
  struct State {
    vec b1, b2;          // [beta, log(scale)] per component
    double s1, s2;       // sigma / phi / shape
    double theta;
    vec gam;
    ivec z;
    vec start1, start2, startg;   // Newton starting points (previous modes)
    double lp;           // log posterior (marginal over z) at the last z-update
  };
  State st;
  st.b1.zeros(D); st.b2.zeros(D);
  st.b1.head(K) = pb1.mean; st.b2.head(K) = pb2.mean;
  st.s1 = st.s2 = 1.0;
  st.theta = theta_a / (theta_a + theta_b);
  st.gam = use_logistic ? pg.mean : vec();
  st.z.ones(N);
  st.lp = NEG_INF;
  if (m.family == GAUSSIAN) st.s1 = st.s2 = sample_sd(y, y_scale);

  // Slice width of the collapsed theta update on the logit scale: about the
  // posterior SD of logit(theta) near the given theta, clamped. During the
  // pilot chains and the burn-in (whose draws are not stored) it follows the
  // current theta, which speeds up convergence; from the first stored sweep on
  // it is fixed, because a width that depends on the current point would make
  // the stepping-out procedure irreversible.
  const double nn_theta = static_cast<double>(nonsafe.n_elem);
  auto theta_width = [&](double th) {
    return std::min(4.0, std::max(0.05, 4.0 / std::sqrt(std::max(1.0, nn_theta * th * (1.0 - th)))));
  };
  double theta_w = theta_width(st.theta);
  bool adapt_width = true;

  // scratch
  vec eta1(N), eta2(N), thn(N);
  Model zmod;   // logistic regression model for gamma | z (Path B)
  mat Zn;       // rows of Z for the non-safe records
  if (use_logistic) {
    zmod.family = BINOMIAL; zmod.X = Z; zmod.N = N; zmod.K = M; zmod.y.set_size(N);
    Zn = Z.rows(nonsafe);
  }
  long acc1 = 0, acc2 = 0, accg = 0, n_mh = 0;
  vec ll1(N), ll2(N);   // per-record component log-likelihoods at the current parameters
  vec dll(nonsafe.n_elem);   // ll1 - ll2 over the non-safe records (collapsed theta update)

  // ---- conditional updates -------------------------------------------------
  // Component log-likelihoods of every record at the current parameter values
  // (shared by the collapsed theta update and the z update).
  auto compute_loglik = [&](const State& s) {
    eta1 = X * s.b1.head(K); eta2 = X * s.b2.head(K);
    if (weibull) { eta1 += s.b1[K]; eta2 += s.b2[K]; }
    for (int n = 0; n < N; ++n) {
      ll1[n] = m.loglik_one(n, eta1[n], s.s1);
      ll2[n] = m.loglik_one(n, eta2[n], s.s2);
      // a NaN (numerical failure) must never decide an allocation: treat it as
      // an impossible observation under that component
      if (std::isnan(ll1[n])) ll1[n] = NEG_INF;
      if (std::isnan(ll2[n])) ll2[n] = NEG_INF;
    }
  };

  // z | rest; also accumulates the log posterior density (likelihood marginal
  // over z plus log priors, up to a constant) of the current parameter values
  // on the natural scale of every parameter. Requires compute_loglik().
  auto update_z = [&](State& s) {
    if (use_logistic) thn = Z * s.gam;
    double lp = 0.0;
    for (int n = 0; n < N; ++n) {
      if (safe[n] != 0) { s.z[n] = 1; lp += ll1[n]; continue; }
      double lt = use_logistic ? -log1p_exp(-thn[n]) : std::log(s.theta);
      double l1t = use_logistic ? -log1p_exp(thn[n]) : std::log1p(-s.theta);
      double a = lt + ll1[n];
      double b = l1t + ll2[n];
      double p1;
      if (!std::isfinite(a) && !std::isfinite(b)) { p1 = 0.5; lp = NEG_INF; }
      else { p1 = inv_logit(a - b); lp += std::max(a, b) + log1p_exp(-std::fabs(a - b)); }
      s.z[n] = (unif_rand() < p1) ? 1 : 2;
    }
    // log priors
    for (int i = 0; i < K; ++i) {
      double u1 = (s.b1[i] - pb1.mean[i]) / pb1.sd[i], u2 = (s.b2[i] - pb2.mean[i]) / pb2.sd[i];
      lp -= 0.5 * (u1 * u1 + u2 * u2);
    }
    if (has_disp) lp += sp1.logd(s.s1) + sp2.logd(s.s2);
    if (weibull) lp += sc1.logd(std::exp(s.b1[K])) + sc2.logd(std::exp(s.b2[K]));
    if (use_logistic) {
      for (int i = 0; i < M; ++i) { double u = (s.gam[i] - pg.mean[i]) / pg.sd[i]; lp -= 0.5 * u * u; }
    } else {
      lp += (theta_a - 1.0) * std::log(s.theta) + (theta_b - 1.0) * std::log1p(-s.theta);
    }
    s.lp = lp;
  };

  // Mixing weight.
  // Path A: theta | beta, dispersions with z integrated out ("collapsed"
  // update) by slice sampling on the logit scale; together with the z draw
  // that follows it this is a blocked update of (theta, z), which mixes far
  // better than the Beta-conjugate draw given z. Requires compute_loglik().
  // Path B: gamma | z by Laplace independence MH (safe records excluded).
  auto update_theta = [&](State& s) {
    if (use_logistic) {
      for (arma::uword i = 0; i < nonsafe.n_elem; ++i) zmod.y[nonsafe[i]] = (s.z[nonsafe[i]] == 1) ? 1.0 : 0.0;
      BlockTarget tg{zmod, Zn, nonsafe, 0.0, pg, vec(), vec()};
      if (laplace_mh(tg, s.gam, s.startg)) ++accg;
    } else if (collapse_theta) {
      // With u = logit(theta): log(theta) = -log1p_exp(-u), log(1 - theta) = -log1p_exp(u),
      // and a record with finite ll1, ll2 contributes log(theta f1 + (1 - theta) f2)
      //   = log(1 - theta) + ll2 + log1p_exp(u + ll1 - ll2); the ll2 term is constant in u.
      // A record that is impossible under component 2 (ll2 = -Inf) contributes
      // log(theta) + ll1, one impossible under component 1 contributes
      // log(1 - theta) + ll2, and one impossible under both carries no
      // information on theta.
      int nf = 0, n_pos = 0, n_neg = 0;
      for (arma::uword i = 0; i < nonsafe.n_elem; ++i) {
        int n = nonsafe[i];
        double d = ll1[n] - ll2[n];
        if (std::isnan(d)) continue;
        if (d == R_PosInf) ++n_pos;
        else if (d == R_NegInf) ++n_neg;
        else dll[nf++] = d;
      }
      const double a_eff = theta_a + n_pos, b_eff = theta_b + nf + n_neg;
      auto logf = [&](double u) {
        double lt = -log1p_exp(-u), l1t = -log1p_exp(u);
        double lp = a_eff * lt + b_eff * l1t;   // Beta(a, b) prior + logit Jacobian
        for (int i = 0; i < nf; ++i) lp += log1p_exp(u + dll[i]);
        return lp;
      };
      // Width: adaptive during pilots and burn-in, fixed while draws are stored
      // (see theta_width above); stepping out adapts to the actual slice.
      double u = slice_sample(std::log(s.theta / (1.0 - s.theta)), logf,
                              adapt_width ? theta_width(s.theta) : theta_w);
      s.theta = std::min(std::max(inv_logit(u), 1e-12), 1.0 - 1e-12);
    } else {
      int n1 = 0, n2 = 0;
      for (arma::uword i = 0; i < nonsafe.n_elem; ++i) { if (s.z[nonsafe[i]] == 1) ++n1; else ++n2; }
      s.theta = R::rbeta(theta_a + n1, theta_b + n2);
      if (s.theta <= 0.0) s.theta = 1e-12;
      if (s.theta >= 1.0) s.theta = 1.0 - 1e-12;
    }
  };

  // component-specific parameters | z
  auto update_components = [&](State& s) {
    for (int j = 1; j <= 2; ++j) {
      vec& b = (j == 1) ? s.b1 : s.b2;
      vec& start = (j == 1) ? s.start1 : s.start2;
      double& sc = (j == 1) ? s.s1 : s.s2;
      const BlockPrior& pb = (j == 1) ? pb1 : pb2;
      const ScalarPrior& sp = (j == 1) ? sp1 : sp2;
      long& acc = (j == 1) ? acc1 : acc2;

      uvec idx = arma::find(s.z == j);
      const int nj = idx.n_elem;
      mat A = X.rows(idx);
      if (weibull) A = arma::join_rows(A, vec(nj, arma::fill::ones));

      if (m.family == GAUSSIAN) {
        gaussian_beta_draw(m, idx, sc, pb, b);
        double rss = 0.0;
        if (nj > 0) { vec r = m.y.elem(idx) - A * b; rss = arma::dot(r, r); }
        auto logf = [&](double u) { return -nj * u - 0.5 * rss * std::exp(-2.0 * u) + sp.logd_log(u); };
        sc = std::exp(slice_sample(std::log(sc), logf));
        continue;
      }

      BlockTarget tb{m, A, idx, sc, pb, vec(), vec()};
      if (laplace_mh(tb, b, start)) ++acc;

      if (m.family == GAMMA) {
        // phi | beta, z through sufficient statistics
        double S1 = 0.0, S2 = 0.0;
        if (nj > 0) {
          vec eta = A * b;
          for (int i = 0; i < nj; ++i) {
            double ly = m.logy[idx[i]];
            S1 += ly - eta[i] - m.y[idx[i]] * std::exp(-eta[i]);
            S2 += ly;
          }
        }
        auto logf = [&](double u) {
          double phi = std::exp(u);
          return nj * (phi * std::log(phi) - std::lgamma(phi)) + phi * S1 - S2 + sp.logd_log(u);
        };
        sc = std::exp(slice_sample(std::log(sc), logf));
      } else if (m.family == SURV_GAMMA) {
        // phi | beta, z: closed form for the events, pgamma for the censored records
        vec eta = A * b;
        double ne = 0.0, S1 = 0.0, S2 = 0.0;
        std::vector<int> cens; std::vector<double> cens_eta;
        for (int i = 0; i < nj; ++i) {
          int n = idx[i];
          if (m.event[n] == 1) { ne += 1.0; S1 += m.logy[n] - eta[i] - m.y[n] * std::exp(-eta[i]); S2 += m.logy[n]; }
          else { cens.push_back(n); cens_eta.push_back(eta[i]); }
        }
        auto logf = [&](double u) {
          double phi = std::exp(u);
          double lp = ne * (phi * std::log(phi) - std::lgamma(phi)) + phi * S1 - S2;
          for (size_t i = 0; i < cens.size(); ++i) lp += gamma_log_surv(phi, m.logy[cens[i]], cens_eta[i]);
          return lp + sp.logd_log(u);
        };
        sc = std::exp(slice_sample(std::log(sc), logf));
      } else if (survival) {
        vec eta = A * b;
        auto logf = [&](double u) {
          double v = std::exp(u);
          return component_loglik(m, idx, eta, v) + sp.logd_log(u);
        };
        sc = std::exp(slice_sample(std::log(sc), logf));
      }
    }
  };

  // ---- collapsed joint moves (partially collapsed Gibbs) ----------------------
  // With few informative observations per record (e.g. binary outcomes) the
  // parameters and the allocations z are strongly dependent, so the updates
  // given z move slowly. All continuous parameters therefore also receive
  // Metropolis moves that target their joint posterior with z integrated out,
  //   phi = (b1, b2, logit(theta) or gamma, log disp1, log disp2),
  //   pi(phi | y) propto prior(phi) * prod_n [theta_n f1(y_n) + (1 - theta_n) f2(y_n)]
  // over the records that are not safe (safe records contribute f1 only;
  // Jacobians of the logit / log transformations included). These moves come
  // after the mixing-weight update and before z is redrawn, so together with
  // the z draw they form a blocked update of (phi, z) (partially collapsed Gibbs
  // sampler; van Dyk and Park, 2008). Two proposals are used: a Gaussian random
  // walk with covariance (2.38^2 / P) C and a multivariate t (8 df) independence
  // proposal with location mu and scale 1.5 C, where mu and C are the mean and
  // covariance of phi over the burn-in draws after the first quarter of the
  // burn-in (adaptive Metropolis, refreshed every 50 draws). From the first
  // stored sweep on they are fixed, so the stored draws come from a fixed, valid
  // kernel. With fewer than max(50, 10 P) such burn-in draws the moves stay off.
  const bool joint_disp = has_disp;
  const int P = 2 * D + (use_logistic ? M : 1) + (joint_disp ? 2 : 0);
  struct JointMove {
    bool active = false, frozen = false;
    long n = 0;                    // burn-in draws collected
    vec mean; mat m2;              // running mean and sum of squared deviations
    vec mu; mat L;                 // proposal location and lower Cholesky factor of C
    long n_rw = 0, acc_rw = 0, n_ind = 0, acc_ind = 0;
  };
  JointMove jm;
  const long jm_min = std::max(50L, 10L * static_cast<long>(P));
  const double JM_IND_SCALE = 1.5;          // inflation of the independence proposal
  vec ll1p(N), ll2p(N), thn_p(N);

  auto pack = [&](const State& s) {
    vec phi(P);
    int k = 0;
    for (int i = 0; i < D; ++i) phi[k++] = s.b1[i];
    for (int i = 0; i < D; ++i) phi[k++] = s.b2[i];
    if (use_logistic) { for (int i = 0; i < M; ++i) phi[k++] = s.gam[i]; }
    else phi[k++] = std::log(s.theta / (1.0 - s.theta));
    if (joint_disp) { phi[k++] = std::log(s.s1); phi[k++] = std::log(s.s2); }
    return phi;
  };
  auto unpack = [&](const vec& phi, State& s) {
    int k = 0;
    for (int i = 0; i < D; ++i) s.b1[i] = phi[k++];
    for (int i = 0; i < D; ++i) s.b2[i] = phi[k++];
    if (use_logistic) { for (int i = 0; i < M; ++i) s.gam[i] = phi[k++]; }
    else s.theta = std::min(std::max(inv_logit(phi[k++]), 1e-12), 1.0 - 1e-12);
    if (joint_disp) { s.s1 = std::exp(phi[k++]); s.s2 = std::exp(phi[k++]); }
  };
  // per-record log-likelihoods of both components at phi
  auto joint_ll = [&](const vec& phi, vec& l1, vec& l2) {
    vec b1 = phi.subvec(0, D - 1), b2 = phi.subvec(D, 2 * D - 1);
    const double d1 = joint_disp ? std::exp(phi[P - 2]) : 1.0;
    const double d2 = joint_disp ? std::exp(phi[P - 1]) : 1.0;
    vec e1 = X * b1.head(K), e2 = X * b2.head(K);
    if (weibull) { e1 += b1[K]; e2 += b2[K]; }
    for (int n = 0; n < N; ++n) {
      double v1 = m.loglik_one(n, e1[n], d1), v2 = m.loglik_one(n, e2[n], d2);
      l1[n] = std::isnan(v1) ? NEG_INF : v1;
      l2[n] = std::isnan(v2) ? NEG_INF : v2;
    }
  };
  // log of the collapsed joint target at phi (up to a constant), given the
  // per-record log-likelihoods of both components at phi
  auto joint_lp = [&](const vec& phi, const vec& l1, const vec& l2) {
    double lp = 0.0;
    for (int i = 0; i < K; ++i) {
      double u1 = (phi[i] - pb1.mean[i]) / pb1.sd[i], u2 = (phi[D + i] - pb2.mean[i]) / pb2.sd[i];
      lp -= 0.5 * (u1 * u1 + u2 * u2);
    }
    if (weibull) lp += sc1.logd_log(phi[K]) + sc2.logd_log(phi[D + K]);
    const int k0 = 2 * D;
    double lt = 0.0, l1t = 0.0;
    if (use_logistic) {
      vec g = phi.subvec(k0, k0 + M - 1);
      for (int i = 0; i < M; ++i) { double u = (g[i] - pg.mean[i]) / pg.sd[i]; lp -= 0.5 * u * u; }
      thn_p = Z * g;
    } else {
      lt = -log1p_exp(-phi[k0]); l1t = -log1p_exp(phi[k0]);
      lp += theta_a * lt + theta_b * l1t;             // Beta prior + logit Jacobian
    }
    if (joint_disp) lp += sp1.logd_log(phi[P - 2]) + sp2.logd_log(phi[P - 1]);
    for (int n = 0; n < N; ++n) {
      if (safe[n] != 0) { lp += l1[n]; continue; }
      double a, c;
      if (use_logistic) { a = -log1p_exp(-thn_p[n]) + l1[n]; c = -log1p_exp(thn_p[n]) + l2[n]; }
      else { a = lt + l1[n]; c = l1t + l2[n]; }
      if (!std::isfinite(a) && !std::isfinite(c)) return NEG_INF;
      lp += std::max(a, c) + log1p_exp(-std::fabs(a - c));
    }
    return std::isnan(lp) ? NEG_INF : lp;
  };
  auto jm_collect = [&](const State& s) {
    vec phi = pack(s);
    if (jm.n == 0) { jm.mean.zeros(P); jm.m2.zeros(P, P); }
    ++jm.n;
    vec delta = phi - jm.mean;
    jm.mean += delta / static_cast<double>(jm.n);
    jm.m2 += delta * (phi - jm.mean).t();
  };
  auto jm_refresh = [&]() {
    if (jm.n < jm_min) return;
    mat C = jm.m2 / static_cast<double>(jm.n - 1);
    mat R;
    if (!C.is_finite() || !chol_pd(C, R)) return;
    jm.mu = jm.mean; jm.L = R.t(); jm.active = true;
  };
  auto joint_move = [&](State& s) {
    if (!jm.active) return;
    vec phi = pack(s);
    double cur = joint_lp(phi, ll1, ll2);             // ll1, ll2 are at the current values
    if (!std::isfinite(cur)) return;
    const double nu = PROPOSAL_DF;
    vec eps(P);
    auto accept = [&](const vec& prop) {
      unpack(prop, s); ll1 = ll1p; ll2 = ll2p; phi = prop;
    };
    // (i) random walk
    for (int i = 0; i < P; ++i) eps[i] = norm_rand();
    vec prop = phi + (2.38 / std::sqrt(static_cast<double>(P))) * (jm.L * eps);
    joint_ll(prop, ll1p, ll2p);
    double pl = joint_lp(prop, ll1p, ll2p);
    ++jm.n_rw;
    if (std::isfinite(pl) && std::log(unif_rand()) < pl - cur) { accept(prop); cur = pl; ++jm.acc_rw; }
    // (ii) independence proposal
    auto logq = [&](const vec& x) {
      vec u = solve_lower(jm.L, x - jm.mu);
      return -0.5 * (nu + P) * std::log1p(arma::dot(u, u) / (JM_IND_SCALE * nu));
    };
    for (int i = 0; i < P; ++i) eps[i] = norm_rand();
    double w = R::rchisq(nu);
    vec prop2 = jm.mu + (jm.L * eps) * std::sqrt(JM_IND_SCALE * nu / w);
    joint_ll(prop2, ll1p, ll2p);
    double pl2 = joint_lp(prop2, ll1p, ll2p);
    ++jm.n_ind;
    if (std::isfinite(pl2) && std::log(unif_rand()) < (pl2 - logq(prop2)) - (cur - logq(phi))) {
      accept(prop2); ++jm.acc_ind;
    }
  };

  // One sweep: component parameters given z, then the mixing weight (Path A:
  // collapsed over z), then the collapsed joint moves, then z given all
  // new parameter values, so that the log posterior recorded by update_z()
  // refers to the stored parameter values.
  auto sweep = [&](State& s) {
    update_components(s); compute_loglik(s); update_theta(s); joint_move(s); update_z(s);
  };

  // ---- initial values ------------------------------------------------------
  // Component 1: mode of its conditional posterior. Known correct matches
  // (safe records) anchor component 1, so they are used on their own when
  // there are enough of them; otherwise all observations are used.
  {
    const uvec safe_idx = arma::find(safe != 0);
    const uvec& init_idx = (static_cast<int>(safe_idx.n_elem) >= D + 2) ? safe_idx : all_idx;
    mat A_init = X.rows(init_idx);
    if (weibull) A_init = arma::join_rows(A_init, vec(init_idx.n_elem, arma::fill::ones));
    BlockTarget t1{m, A_init, init_idx, st.s1, pb1, vec(), vec()};
    mat P; double lp; vec mode = st.b1;
    if (find_mode(t1, mode, P, lp)) st.b1 = mode;
    if (m.family == GAUSSIAN) st.s1 = sample_sd(y.elem(init_idx) - A_init * st.b1, y_scale);
  }
  // Component 2 (mismatches): the marginal outcome distribution, i.e. an
  // intercept-only model on all observations (slopes at their prior mean).
  // Without an intercept column the coefficients stay at their prior mean
  // (for the Weibull model the log-scale takes the marginal value instead).
  double marginal_value = has_intercept ? pb2.mean[0] : 0.0;
  bool have_marginal = false;
  {
    double ybar = arma::mean(y);
    switch (m.family) {
      case GAUSSIAN: marginal_value = ybar; have_marginal = true; break;
      case BINOMIAL: {
        double p = std::min(std::max(ybar, 0.02), 0.98);
        marginal_value = std::log(p / (1.0 - p)); have_marginal = true; break;
      }
      default:                                        // log link / AFT
        if (ybar > 0.0 && std::isfinite(std::log(ybar))) { marginal_value = std::log(ybar); have_marginal = true; }
        break;
    }
    // respect a tight prior on the intercept
    if (have_marginal && has_intercept) {
      double u = (marginal_value - pb2.mean[0]) / pb2.sd[0];
      if (std::fabs(u) > 3.0) marginal_value = pb2.mean[0] + 3.0 * (u > 0 ? 1.0 : -1.0) * pb2.sd[0];
    }
  }
  // starting values of the mismatch component at the marginal distribution
  auto set_mismatch_start = [&](State& s) {
    s.b2.head(K) = pb2.mean;
    if (have_marginal) {
      if (has_intercept) s.b2[0] = marginal_value;
      else if (weibull) s.b2[K] = marginal_value;
    }
    if (m.family == GAUSSIAN) s.s2 = y_scale;
    s.start2 = s.b2;
  };

  const bool user_init = init.length() > 0;
  if (init.containsElementNamed("beta1"))  st.b1.head(K) = as<vec>(init["beta1"]);
  if (init.containsElementNamed("beta2"))  st.b2.head(K) = as<vec>(init["beta2"]);
  if (init.containsElementNamed("disp1"))  st.s1 = as<double>(init["disp1"]);
  if (init.containsElementNamed("disp2"))  st.s2 = as<double>(init["disp2"]);
  if (weibull && init.containsElementNamed("scale1")) st.b1[K] = std::log(as<double>(init["scale1"]));
  if (weibull && init.containsElementNamed("scale2")) st.b2[K] = std::log(as<double>(init["scale2"]));
  if (!use_logistic && init.containsElementNamed("theta")) st.theta = as<double>(init["theta"]);
  if (use_logistic && init.containsElementNamed("gamma")) st.gam = as<vec>(init["gamma"]);
  if (!(st.s1 > 0.0) || !(st.s2 > 0.0)) stop("Initial dispersion parameters must be positive.");
  if (!(st.theta > 0.0 && st.theta < 1.0)) stop("The initial value of theta must lie in (0, 1).");
  if (!st.b1.is_finite() || !st.b2.is_finite() || (use_logistic && !st.gam.is_finite()))
    stop("Initial values must be finite.");
  st.start1 = st.b1; st.start2 = st.b2; st.startg = st.gam;

  // Interrupts (Esc / Ctrl-C) are honoured at most every quarter of a second,
  // whatever the cost of a sweep.
  auto last_check = std::chrono::steady_clock::now();
  auto check_interrupt = [&]() {
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration<double>(now - last_check).count() > 0.25) {
      Rcpp::checkUserInterrupt();
      last_check = now;
    }
  };

  // Change in the log posterior (with z integrated out) when the two labels
  // of a state are exchanged (b1 <-> b2 including the Weibull log-scale,
  // s1 <-> s2, theta -> 1 - theta or gamma -> -gamma). It is used only
  // without safe matches, so the likelihood with z integrated out is unchanged
  // and only the priors contribute (on any scale of the parameters: the
  // Jacobians of the log and logit transformations are symmetric under the
  // exchange): the component-specific priors and, when it is not symmetric
  // (pre_orient = 2), the prior on the match probability. Computed term by
  // term, so that it is exactly 0 when the component-specific priors are
  // identical and the prior on the match probability is symmetric.
  auto exchange_lp_change = [&](const State& s) {
    double d = 0.0;
    for (int i = 0; i < K; ++i) {
      double a1 = (s.b2[i] - pb1.mean[i]) / pb1.sd[i], a2 = (s.b1[i] - pb2.mean[i]) / pb2.sd[i];
      double c1 = (s.b1[i] - pb1.mean[i]) / pb1.sd[i], c2 = (s.b2[i] - pb2.mean[i]) / pb2.sd[i];
      d -= 0.5 * ((a1 * a1 + a2 * a2) - (c1 * c1 + c2 * c2));
    }
    if (has_disp) d += (sp1.logd(s.s2) + sp2.logd(s.s1)) - (sp1.logd(s.s1) + sp2.logd(s.s2));
    if (weibull) {
      const double e1 = std::exp(s.b1[K]), e2 = std::exp(s.b2[K]);
      d += (sc1.logd(e2) + sc2.logd(e1)) - (sc1.logd(e1) + sc2.logd(e2));
    }
    if (use_logistic) {
      for (int i = 0; i < M; ++i) {
        double u = (-s.gam[i] - pg.mean[i]) / pg.sd[i], v = (s.gam[i] - pg.mean[i]) / pg.sd[i];
        d -= 0.5 * (u * u - v * v);
      }
    } else if (theta_a != theta_b) {
      d += (theta_a - theta_b) * (std::log1p(-s.theta) - std::log(s.theta));
    }
    return d;
  };
  // whether the orientation of the starting state can apply, and by which rule
  const bool orient_possible = pre_orient != 0 && !user_init && nonsafe.n_elem == static_cast<arma::uword>(N);
  const bool orient_by_mass = pre_orient == 2;
  // share of the records not flagged as safe allocated to component 1 in a state
  auto share1_of = [&](const State& s) {
    arma::uword c1 = 0;
    for (arma::uword i = 0; i < nonsafe.n_elem; ++i) if (s.z[nonsafe[i]] == 1) ++c1;
    return nonsafe.n_elem > 0 ? static_cast<double>(c1) / static_cast<double>(nonsafe.n_elem) : NA_REAL;
  };

  // ---- multi-start pilot runs ---------------------------------------------
  // Mixture posteriors are multimodal and a Gibbs sampler cannot jump between
  // modes. Several short pilot chains are therefore run from different
  // starting configurations and the main chain continues from the one with
  // the highest average log posterior. Each pilot chain also records, over
  // the same final sweeps, the average share of the records not flagged as
  // safe allocated to component 1 (which decides the orientation under the
  // majority convention: the allocation of a single sweep can fall on the
  // wrong side of one half when the components are nearly balanced) and,
  // when the starting state may be oriented, the log of the mean of
  // exp(exchange_lp_change): the log ratio of the posterior probabilities of
  // the exchanged labelling and of the labelling of the chain, since the
  // posterior probability of the exchanged region is the posterior
  // expectation of exp(exchange_lp_change) over the region of the chain.
  NumericVector pilot_lp(n_pilot > 0 && !user_init ? n_pilot : 0);
  double best_exchange = NA_REAL;   // estimate of the selected pilot chain
  double best_share1 = NA_REAL;     // average share of component 1 in it
  if (n_pilot > 0 && !user_init) {
    const int n_avg = std::max(1, std::min(50, pilot_iter / 2));
    State best = st; double best_score = NEG_INF;
    for (int p = 0; p < n_pilot; ++p) {
      State s = st;
      if (p <= 1) {
        if (p == 1) set_mismatch_start(s);   // mismatch component = marginal distribution
        compute_loglik(s); update_z(s);      // allocations from the initial parameters
      } else {                               // random allocation
        set_mismatch_start(s);
        if (use_logistic) thn = Z * s.gam;
        for (int n = 0; n < N; ++n) {
          double pr = use_logistic ? inv_logit(thn[n]) : s.theta;
          s.z[n] = (safe[n] != 0 || unif_rand() < pr) ? 1 : 2;
        }
      }
      double score = 0.0, share = 0.0;
      LogMeanExp exch;
      for (int it = 1; it <= pilot_iter; ++it) {
        check_interrupt();
        sweep(s);
        if (it > pilot_iter - n_avg) {
          score += s.lp / n_avg;
          if (nonsafe.n_elem > 0) share += share1_of(s) / n_avg;
          if (orient_possible) exch.add(exchange_lp_change(s));
        }
      }
      pilot_lp[p] = score;
      if (verbose) Rprintf("Pilot chain %d / %d: average log posterior %.2f\n", p + 1, n_pilot, score);
      if (std::isfinite(score) && score > best_score) {
        best_score = score; best = s;
        best_exchange = exch.value();
        best_share1 = (nonsafe.n_elem > 0) ? share : NA_REAL;
      }
    }
    if (std::isfinite(best_score)) {
      st = best;
    } else {
      best_exchange = NA_REAL;
      best_share1 = NA_REAL;
      // no pilot reached a finite log posterior (reported by the R code):
      // start the main chain from the initial values
      compute_loglik(st); update_z(st);
    }
    acc1 = acc2 = accg = 0;
  } else {
    compute_loglik(st); update_z(st);
  }

  // ---- orientation of the starting state ----------------------------------
  // The two labels of the starting state of the main chain (the selected
  // pilot state, or the initial allocation) can be exchanged: the coefficient
  // blocks (including the Weibull log-scale), the dispersions, theta -> 1 -
  // theta or gamma -> -gamma, the allocations and every cached quantity
  // (Newton starting points, and the component log-likelihoods and linear
  // predictors, which the next sweep recomputes before using them). This
  // changes only the state the main chain starts from, so the transition
  // kernel of the main chain and its stationary distribution are unchanged.
  // The collapsed joint moves have collected nothing yet (they adapt during
  // the main burn-in), and the log posterior of the state is recomputed by the
  // next sweep. User starting values are respected, and safe matches (fixed in
  // component 1) identify the labels, so neither is ever reoriented.
  //
  // pre_orient = 2 (the prior on the match probability identifies the
  // labels): the labelling with the larger posterior probability is wanted.
  // A Gibbs chain stays in the labelling it starts in, which the pilot chains
  // reach more or less at random when the two labellings fit the data equally
  // well, so the labels are exchanged when the estimated log ratio of the
  // posterior probabilities of the exchanged and the original labelling
  // (from the selected pilot chain; exchange_lp_change() at the starting
  // state when that estimate is missing or NaN, e.g. without pilot chains) is
  // positive. A chain started in the mirror image of a region where a pilot
  // chain stayed is equally stable there, so no check chains are needed.
  //
  // pre_orient = 1 (majority convention: neither safe matches nor the prior
  // on the match probability identify the labels): component 1 is meant to be
  // the majority component. The labels are exchanged when fewer than half of
  // the records are allocated to component 1 in the starting state, on
  // average over the last sweeps of the selected pilot chain (the allocation
  // of a single sweep can fall on the wrong side of one half when the
  // components are nearly balanced; without pilot chains the initial
  // allocation decides, which only approximates the labelling that the chain
  // will reach).
  // Under component-specific priors that differ, the exchanged labelling can
  // be far less probable than the labelling of the starting state (e.g. a
  // tight prior on the slopes of component 2 that the majority of the records
  // contradict), and a Gibbs chain started there may not return within the
  // run. The estimated log ratio of the posterior probabilities of the two
  // labellings (from the selected pilot chain; exchange_lp_change() at the
  // starting state when that estimate is missing or NaN, e.g. without pilot
  // chains) is exactly 0 for identical component-specific priors, and the
  // labels are exchanged when it is at least -orient_tol (a NaN estimate
  // leaves the majority convention in force). Below -orient_tol the estimate
  // alone does not decide: it refers to the mirror image of the region that
  // the pilot chain visited, and a chain started there can move on to another
  // mode of the same labelling that is as probable as the original one (e.g.
  // logistic components with many covariates, whose slopes the default priors
  // shrink differently in the two components). With pilot chains, two check
  // chains of pilot_iter sweeps are therefore run from the selected pilot
  // state, one with the labels exchanged and one without, and the exchange is
  // refused (`orient_refused`: the component-specific priors identify the
  // labels) when, over the second half of its sweeps, the exchanged chain
  // returned to the original labelling (component 1 holds fewer than half of
  // the records on average), or its average log posterior stays below that of
  // the other chain by more than the margin max(orient_tol, 2 se), where se is
  // the Monte Carlo standard error of the difference of the two averages:
  // se^2 = se_1^2 + se_2^2 with se_j^2 = max(V - v_j, v_j / n) for an average
  // over n consecutive sweeps whose log posterior has variance v_j, V being
  // the variance of the log posterior over the posterior. By the law of total
  // variance, the variance of such an average is V minus the expected
  // variance within the n sweeps: a slowly mixing chain explores only part of
  // the posterior within them (v_j < V), so its average is as uncertain as
  // the unexplored part of V, while v_j / n is the error of a chain that mixes
  // well (or drifts, v_j > V). V is the larger of P / 2, its value under a
  // normal approximation of the posterior of the P continuous parameters, and
  // the variance of the log posterior over all sweeps of the chain started
  // without exchange. A level difference is thus judged against the tolerance
  // itself when the chains mix well, and against a few standard deviations of
  // the log posterior when they mix slowly (logistic components with many
  // covariates).
  // Without pilot chains, or for an estimate of -Inf (an exchanged labelling
  // of probability 0), the estimate decides alone. The check chains only
  // choose the state the main chain starts from (they store nothing).
  auto exchange_labels = [&](State& s) {
    s.b1.swap(s.b2);
    s.start1.swap(s.start2);
    std::swap(s.s1, s.s2);
    if (use_logistic) {
      s.gam = -s.gam;
      s.startg = -s.startg;
    } else {
      s.theta = 1.0 - s.theta;
    }
    for (int n = 0; n < N; ++n) s.z[n] = 3 - s.z[n];
  };
  // a check chain from `s`: average and variance of the log posterior, and
  // average share of the records not flagged as safe allocated to component
  // 1, over the last n_chk of pilot_iter sweeps (the average log posterior is
  // not finite when a log posterior is not), and the variance of the log
  // posterior over all pilot_iter sweeps
  const int n_chk = std::max(1, pilot_iter / 2);
  struct Running {                     // running mean and variance (Welford)
    double mean = 0.0, m2 = 0.0; long k = 0;
    void add(double x) {
      if (!std::isfinite(x)) return;
      ++k; const double d = x - mean; mean += d / static_cast<double>(k); m2 += d * (x - mean);
    }
    double var() const { return (k > 1) ? m2 / static_cast<double>(k - 1) : 0.0; }
  };
  auto check_chain = [&](State s, double& avg, double& var, double& share1, double& var_all) {
    Running last, all;
    double sum = 0.0, sh = 0.0;
    for (int it = 1; it <= pilot_iter; ++it) {
      check_interrupt();
      sweep(s);
      all.add(s.lp);
      if (it > pilot_iter - n_chk) {
        sum += s.lp / n_chk;
        last.add(s.lp);
        sh += share1_of(s) / n_chk;
      }
    }
    avg = sum;
    var = last.var();
    share1 = sh;
    var_all = all.var();
  };
  bool pre_oriented = false, orient_refused = false;
  double start_share1 = NA_REAL, orient_change = NA_REAL;
  double check_exch = NA_REAL, check_kept = NA_REAL, check_margin = NA_REAL, check_share1 = NA_REAL;
  if (nonsafe.n_elem > 0) {
    // share of component 1 in the starting state: on average over the last
    // sweeps of the selected pilot chain, else in the initial allocation
    start_share1 = std::isnan(best_share1) ? share1_of(st) : best_share1;
    bool exchange = false;
    if (orient_possible && orient_by_mass) {
      orient_change = std::isnan(best_exchange) ? exchange_lp_change(st) : best_exchange;
      exchange = orient_change > 0.0;          // NaN compares false: no exchange
    } else if (orient_possible && start_share1 < 0.5) {
      orient_change = std::isnan(best_exchange) ? exchange_lp_change(st) : best_exchange;
      // NaN (no estimate) compares false: the majority convention stays in force
      orient_refused = orient_change < -orient_tol;
      if (orient_refused && n_pilot > 0 && std::isfinite(orient_change)) {
        State ex = st;
        exchange_labels(ex);
        double v_exch = 0.0, v_kept = 0.0, share_kept = 0.0, va_exch = 0.0, va_kept = 0.0;
        check_chain(ex, check_exch, v_exch, check_share1, va_exch);
        check_chain(st, check_kept, v_kept, share_kept, va_kept);
        const double V = std::max(0.5 * static_cast<double>(P), va_kept);
        auto se2 = [&](double v) { return std::max(V - v, v / static_cast<double>(n_chk)); };
        check_margin = std::max(orient_tol, 2.0 * std::sqrt(se2(v_exch) + se2(v_kept)));
        // refused when the exchanged chain returned to the original labelling,
        // or when the kept chain has a finite level that the exchanged chain
        // stays clearly below (a non-finite exchanged level counts as below)
        orient_refused = 2.0 * check_share1 < 1.0 ||
          (std::isfinite(check_kept) && !(check_exch >= check_kept - check_margin));
        acc1 = acc2 = accg = 0;
        // the caches refer to the last check sweep: recompute them for `st`
        compute_loglik(st);
        if (use_logistic) thn = Z * st.gam;
        if (verbose)
          Rprintf("Orientation check: average log posterior %.2f with the labels exchanged (share of "
                  "component 1 %.2f), %.2f without (margin %.2f).\n", check_exch, check_share1,
                  check_kept, check_margin);
      }
      exchange = !orient_refused;
    }
    if (exchange) {
      exchange_labels(st);
      if (use_logistic) thn = -thn;
      ll1.swap(ll2);
      eta1.swap(eta2);
      pre_oriented = true;
      if (verbose) {
        if (orient_by_mass)
          Rprintf("Starting state oriented: the labelling with the components exchanged is more probable "
                  "(log ratio %.2f).\n", orient_change);
        else
          Rprintf("Starting state oriented: component 1 now holds the majority of the records.\n");
      }
    }
    if (orient_refused && verbose)
      Rprintf("Starting state not oriented: the component-specific priors make the labelling with "
              "component 1 as the majority less probable (log ratio %.1f).\n", orient_change);
  }

  // ---- storage -------------------------------------------------------------
  // exch_s: exchange_lp_change() of every stored draw, the change in its log
  // posterior when its labels are exchanged (NA with safe matches, for which
  // the likelihood changes too); the R code uses it for the log posterior of
  // relabelled draws and to estimate the posterior probability of the other
  // labelling
  mat beta1_s(S, K), beta2_s(S, K), gamma_s(S, M);
  vec theta_s(S), disp1_s(S), disp2_s(S), scale1_s(S), scale2_s(S), lp_s(S);
  vec exch_s(S);
  exch_s.fill(NA_REAL);
  const bool store_exchange = nonsafe.n_elem == static_cast<arma::uword>(N);

  int stored = 0;
  const int report = std::max(1, n_iter / 10);

  // ---- main chain ----------------------------------------------------------
  for (int it = 1; it <= n_iter; ++it) {
    check_interrupt();
    if (verbose && (it % report == 0 || it == n_iter))
      Rprintf("Iteration %d / %d%s\n", it, n_iter, it <= n_burnin ? " (burn-in)" : "");
    if (it == n_burnin + 1) {          // from the first stored sweep on, fixed proposals
      theta_w = theta_width(st.theta);
      adapt_width = false;
      jm.active = false;
      jm_refresh();                      // final estimate from all burn-in draws (or stays off)
      jm.frozen = jm.active;
      jm.n_rw = jm.acc_rw = jm.n_ind = jm.acc_ind = 0;
      // the acceptance rates of the blocks, like those of the joint moves,
      // refer to the sweeps after the burn-in (all sweeps when n_burnin = 0)
      acc1 = acc2 = accg = 0;
      n_mh = 0;
    }

    sweep(st);
    ++n_mh;

    // adaptation of the collapsed joint moves during the burn-in (after its first quarter)
    if (it <= n_burnin && it > n_burnin / 4) {
      jm_collect(st);
      if (jm.n % 50 == 0) jm_refresh();
    }

    if (it > n_burnin && (it - n_burnin) % thin == 0 && stored < S) {
      beta1_s.row(stored) = st.b1.head(K).t();
      beta2_s.row(stored) = st.b2.head(K).t();
      if (use_logistic) gamma_s.row(stored) = st.gam.t(); else theta_s[stored] = st.theta;
      disp1_s[stored] = st.s1; disp2_s[stored] = st.s2;
      if (weibull) { scale1_s[stored] = std::exp(st.b1[K]); scale2_s[stored] = std::exp(st.b2[K]); }
      lp_s[stored] = st.lp;
      if (store_exchange) exch_s[stored] = exchange_lp_change(st);
      for (int n = 0; n < N; ++n) z_s(stored, n) = st.z[n];
      ++stored;
    }
  }

  List out = List::create(
    _["beta1"] = beta1_s,
    _["beta2"] = beta2_s,
    _["z"] = z_s,
    _["lp"] = lp_s,
    _["exchange_lp"] = exch_s,
    _["n_draws"] = stored
  );
  if (use_logistic) out["gamma"] = gamma_s; else out["theta"] = theta_s;
  if (has_disp) { out["disp1"] = disp1_s; out["disp2"] = disp2_s; }
  if (weibull) { out["scale1"] = scale1_s; out["scale2"] = scale2_s; }
  // acceptance rates of the blocks over the sweeps after the burn-in (over all
  // sweeps without burn-in), the window of the joint moves below
  out["accept"] = List::create(
    _["beta1"] = (m.family == GAUSSIAN) ? 1.0 : static_cast<double>(acc1) / n_mh,
    _["beta2"] = (m.family == GAUSSIAN) ? 1.0 : static_cast<double>(acc2) / n_mh,
    _["gamma"] = use_logistic ? static_cast<double>(accg) / n_mh : NA_REAL
  );
  out["pilot_lp"] = pilot_lp;
  // orientation before the main chain: whether the labels of the starting
  // state were exchanged, the share of the records not flagged as safe that
  // the starting state allocated to component 1 before that (on average over
  // the last sweeps of the selected pilot chain, else in the initial
  // allocation), whether an exchange was refused because the
  // component-specific priors make the exchanged labelling much less probable
  // (majority convention), the estimated log ratio of the posterior
  // probabilities of the exchanged and the original labelling (NA when no
  // exchange was considered), and the check chains run when that
  // estimate was below -orient_tol (average log posteriors with the labels
  // exchanged and kept, the margin, and the average share of the records
  // allocated to component 1 in the exchanged chain; NA when no check was run)
  out["pre_oriented"] = pre_oriented;
  out["start_share1"] = start_share1;
  out["orient_refused"] = orient_refused;
  out["orient_lp_change"] = orient_change;
  out["orient_check"] = NumericVector::create(_["exchanged"] = check_exch, _["kept"] = check_kept,
                                              _["margin"] = check_margin, _["share1"] = check_share1);
  // acceptance rates of the collapsed moves while draws were stored (NA when off)
  auto rate = [](long a, long n) { return n > 0 ? static_cast<double>(a) / n : NA_REAL; };
  out["accept_collapsed"] = List::create(
    _["random_walk"] = jm.frozen ? rate(jm.acc_rw, jm.n_rw) : NA_REAL,
    _["independence"] = jm.frozen ? rate(jm.acc_ind, jm.n_ind) : NA_REAL
  );
  // burn-in draws collected for the proposal of the joint moves and the number
  // they need (the moves are off with fewer, e.g. after a short burn-in)
  out["joint_burnin"] = NumericVector::create(_["collected"] = static_cast<double>(jm.n),
                                              _["needed"] = static_cast<double>(jm_min));
  return out;
}
