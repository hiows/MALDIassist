#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace spectrum_math {

constexpr double SQRT_2PI = 2.50662827463100050242;
constexpr int DEFAULT_MAX_ITER = 100;
constexpr double DEFAULT_KDE_CUTOFF = 5.0;

inline bool is_finite_double(const double x) {
  return R_finite(x);
}

inline void validate_numeric_vector(const NumericVector& v, const char* name,
                                    const bool allow_empty = false) {
  if (!allow_empty && v.size() == 0) {
    stop(std::string(name) + " must not be empty.");
  }
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    if (!is_finite_double(v[i])) {
      stop(std::string(name) + " contains NA, NaN, Inf, or -Inf.");
    }
  }
}

inline void validate_kde_inputs(const NumericVector& x,
                                const NumericVector& x_obs,
                                const NumericVector& y_obs,
                                const double bw) {
  validate_numeric_vector(x, "x", true);
  validate_numeric_vector(x_obs, "x_obs");
  validate_numeric_vector(y_obs, "y_obs");

  if (x_obs.size() != y_obs.size()) {
    stop("x_obs and y_obs must have the same length.");
  }
  if (!is_finite_double(bw) || bw <= 0.0) {
    stop("bw must be a finite positive number.");
  }
}

inline bool is_strictly_increasing(const NumericVector& v) {
  if (v.size() < 2) return true;
  for (R_xlen_t i = 1; i < v.size(); ++i) {
    if (v[i] <= v[i - 1]) return false;
  }
  return true;
}

inline std::pair<R_xlen_t, R_xlen_t> kde_support_bounds(
    const double xi,
    const NumericVector& x_obs,
    const double radius) {
  const R_xlen_t n_obs = x_obs.size();
  const double lo = xi - radius;
  const double hi = xi + radius;

  R_xlen_t j_lo = 0;
  while (j_lo < n_obs && x_obs[j_lo] < lo) {
    ++j_lo;
  }

  R_xlen_t j_hi = n_obs - 1;
  while (j_hi >= j_lo && x_obs[j_hi] > hi) {
    --j_hi;
  }

  if (j_lo > j_hi) {
    return std::make_pair(n_obs, n_obs - 1);
  }
  return std::make_pair(j_lo, j_hi);
}

inline void kde_accumulate(const double xi,
                           const NumericVector& x_obs,
                           const NumericVector& y_obs,
                           const double bw,
                           const R_xlen_t j_lo,
                           const R_xlen_t j_hi,
                           double& n0,
                           double& d0,
                           double& n1,
                           double& d1,
                           double& n2,
                           double& d2,
                           double& n3,
                           double& d3) {
  const double bw2 = bw * bw;
  const double bw4 = bw2 * bw2;
  const double bw6 = bw4 * bw2;
  const double inv_bw2 = 1.0 / bw2;
  const double inv_bw4 = 1.0 / bw4;
  const double inv_bw6 = 1.0 / bw6;
  const double kernel_const = 1.0 / (bw * SQRT_2PI);

  n0 = 0.0; d0 = 0.0;
  n1 = 0.0; d1 = 0.0;
  n2 = 0.0; d2 = 0.0;
  n3 = 0.0; d3 = 0.0;

  for (R_xlen_t j = j_lo; j <= j_hi; ++j) {
    const double diff = xi - x_obs[j];
    const double yj = y_obs[j];
    const double k = std::exp(-0.5 * diff * diff * inv_bw2) * kernel_const;
    const double k1 = -diff * inv_bw2 * k;
    const double k2 = (diff * diff * inv_bw4 - inv_bw2) * k;
    const double k3 = diff * (3.0 * bw2 - diff * diff) * inv_bw6 * k;

    n0 += k * yj;  d0 += k;
    n1 += k1 * yj; d1 += k1;
    n2 += k2 * yj; d2 += k2;
    n3 += k3 * yj; d3 += k3;
  }
}

inline double kde_deriv_scalar(const double xi,
                               const NumericVector& x_obs,
                               const NumericVector& y_obs,
                               const double bw,
                               const int deriv_order,
                               const bool use_truncation,
                               const double cutoff) {
  const double denom_floor = std::numeric_limits<double>::min();

  R_xlen_t j_lo = 0;
  R_xlen_t j_hi = x_obs.size() - 1;

  if (use_truncation) {
    const std::pair<R_xlen_t, R_xlen_t> bounds =
      kde_support_bounds(xi, x_obs, cutoff * bw);
    j_lo = bounds.first;
    j_hi = bounds.second;
    if (j_lo > j_hi) {
      return NA_REAL;
    }
  }

  double n0 = 0.0, d0 = 0.0;
  double n1 = 0.0, d1 = 0.0;
  double n2 = 0.0, d2 = 0.0;
  double n3 = 0.0, d3 = 0.0;

  kde_accumulate(xi, x_obs, y_obs, bw, j_lo, j_hi,
                 n0, d0, n1, d1, n2, d2, n3, d3);

  if (d0 <= denom_floor) {
    return NA_REAL;
  }

  if (deriv_order == 0) {
    return n0 / d0;
  }
  if (deriv_order == 1) {
    return (n1 * d0 - n0 * d1) / (d0 * d0);
  }
  if (deriv_order == 2) {
    const double d0_2 = d0 * d0;
    const double d0_3 = d0_2 * d0;
    return n2 / d0 - n0 * d2 / d0_2 - 2.0 * n1 * d1 / d0_2 +
      2.0 * n0 * d1 * d1 / d0_3;
  }

  const double d0_2 = d0 * d0;
  const double d0_3 = d0_2 * d0;
  const double d0_4 = d0_3 * d0;
  return n3 / d0 -
    3.0 * n2 * d1 / d0_2 +
    3.0 * n1 * (2.0 * d1 * d1 / d0_3 - d2 / d0_2) +
    n0 * (-d3 / d0_2 + 6.0 * d1 * d2 / d0_3 -
          6.0 * d1 * d1 * d1 / d0_4);
}

inline void validate_root_inputs(const NumericVector& x,
                                 const double yi,
                                 const double tol,
                                 const int max_iter) {
  validate_numeric_vector(x, "x");
  if (x.size() < 2) {
    stop("x must contain at least two points.");
  }
  if (!is_finite_double(yi)) {
    stop("yi must be a finite number.");
  }
  if (!is_finite_double(tol) || tol <= 0.0) {
    stop("tol must be a finite positive number.");
  }
  if (max_iter <= 0) {
    stop("max_iter must be a positive integer.");
  }

  bool increasing = true;
  bool decreasing = true;
  for (R_xlen_t i = 1; i < x.size(); ++i) {
    increasing = increasing && (x[i] > x[i - 1]);
    decreasing = decreasing && (x[i] < x[i - 1]);
  }
  if (!increasing && !decreasing) {
    stop("x must be strictly monotonic for interval-based root finding.");
  }
}

inline double eval_scalar_function(const Function& fun, const double x) {
  NumericVector y = fun(NumericVector::create(x));
  if (y.size() != 1) {
    stop("fun must return a numeric vector with the same length as its input.");
  }
  if (!is_finite_double(y[0])) {
    stop("fun returned NA, NaN, Inf, or -Inf during bisection.");
  }
  return y[0];
}

inline NumericVector eval_vector_function(const Function& fun, const NumericVector& x) {
  NumericVector y = fun(x);
  if (y.size() != x.size()) {
    stop("fun must return a numeric vector with the same length as x.");
  }
  validate_numeric_vector(y, "fun(x)");
  return y;
}

inline bool close_to_existing_root(const std::vector<double>& roots,
                                   const double candidate,
                                   const double tol) {
  for (double root : roots) {
    if (std::abs(root - candidate) <= tol) return true;
  }
  return false;
}

inline double bisection_root_cached(const Function& fun,
                                    double left,
                                    double right,
                                    const double yi,
                                    const double tol,
                                    const int max_iter) {
  double f_left = eval_scalar_function(fun, left) - yi;
  double f_right = eval_scalar_function(fun, right) - yi;

  if (std::abs(f_left) <= tol) return left;
  if (std::abs(f_right) <= tol) return right;
  if (f_left * f_right > 0.0) {
    stop("Bisection requires a sign-changing interval after subtracting yi.");
  }

  double mid = 0.5 * (left + right);
  double f_mid = eval_scalar_function(fun, mid) - yi;

  for (int iter = 0; iter < max_iter; ++iter) {
    mid = 0.5 * (left + right);
    f_mid = eval_scalar_function(fun, mid) - yi;

    if (std::abs(f_mid) <= tol || std::abs(right - left) <= tol) {
      return mid;
    }

    if (f_left * f_mid <= 0.0) {
      right = mid;
      f_right = f_mid;
    } else {
      left = mid;
      f_left = f_mid;
    }
  }

  return mid;
}

inline double bisection_root_kde_deriv(const NumericVector& x_obs,
                                       const NumericVector& y_obs,
                                       const double bw,
                                       const int deriv_order,
                                       const bool use_truncation,
                                       const double cutoff,
                                       double left,
                                       double right,
                                       const double yi,
                                       const double tol,
                                       const int max_iter) {
  double f_left = kde_deriv_scalar(left, x_obs, y_obs, bw, deriv_order,
                                   use_truncation, cutoff) - yi;
  double f_right = kde_deriv_scalar(right, x_obs, y_obs, bw, deriv_order,
                                    use_truncation, cutoff) - yi;

  if (!is_finite_double(f_left) || !is_finite_double(f_right)) {
    stop("KDE derivative evaluation returned non-finite values during bisection.");
  }

  if (std::abs(f_left) <= tol) return left;
  if (std::abs(f_right) <= tol) return right;
  if (f_left * f_right > 0.0) {
    stop("Bisection requires a sign-changing interval after subtracting yi.");
  }

  double mid = 0.5 * (left + right);
  double f_mid = kde_deriv_scalar(mid, x_obs, y_obs, bw, deriv_order,
                                  use_truncation, cutoff) - yi;

  for (int iter = 0; iter < max_iter; ++iter) {
    mid = 0.5 * (left + right);
    f_mid = kde_deriv_scalar(mid, x_obs, y_obs, bw, deriv_order,
                             use_truncation, cutoff) - yi;

    if (!is_finite_double(f_mid)) {
      stop("KDE derivative evaluation returned non-finite values during bisection.");
    }

    if (std::abs(f_mid) <= tol || std::abs(right - left) <= tol) {
      return mid;
    }

    if (f_left * f_mid <= 0.0) {
      right = mid;
      f_right = f_mid;
    } else {
      left = mid;
      f_left = f_mid;
    }
  }

  return mid;
}

inline NumericVector gauss_kde_eval(const NumericVector& x,
                                    const NumericVector& x_obs,
                                    const NumericVector& y_obs,
                                    const double bw,
                                    const int deriv_order,
                                    const bool use_truncation = true,
                                    const double cutoff = DEFAULT_KDE_CUTOFF) {
  validate_kde_inputs(x, x_obs, y_obs, bw);
  if (deriv_order < 0 || deriv_order > 3) {
    stop("deriv_order must be 0, 1, 2, or 3.");
  }

  const bool sorted_obs = is_strictly_increasing(x_obs);
  const bool truncate = use_truncation && sorted_obs;
  const R_xlen_t n = x.size();
  const R_xlen_t n_obs = x_obs.size();
  NumericVector out(n);
  const double denom_floor = std::numeric_limits<double>::min();

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];

    R_xlen_t j_lo = 0;
    R_xlen_t j_hi = n_obs - 1;

    if (truncate) {
      const std::pair<R_xlen_t, R_xlen_t> bounds =
        kde_support_bounds(xi, x_obs, cutoff * bw);
      j_lo = bounds.first;
      j_hi = bounds.second;
      if (j_lo > j_hi) {
        out[i] = NA_REAL;
        continue;
      }
    }

    double n0 = 0.0, d0 = 0.0;
    double n1 = 0.0, d1 = 0.0;
    double n2 = 0.0, d2 = 0.0;
    double n3 = 0.0, d3 = 0.0;

    kde_accumulate(xi, x_obs, y_obs, bw, j_lo, j_hi,
                   n0, d0, n1, d1, n2, d2, n3, d3);

    if (d0 <= denom_floor) {
      out[i] = NA_REAL;
      continue;
    }

    if (deriv_order == 0) {
      out[i] = n0 / d0;
    } else if (deriv_order == 1) {
      out[i] = (n1 * d0 - n0 * d1) / (d0 * d0);
    } else if (deriv_order == 2) {
      const double d0_2 = d0 * d0;
      const double d0_3 = d0_2 * d0;
      out[i] = n2 / d0 - n0 * d2 / d0_2 - 2.0 * n1 * d1 / d0_2 +
        2.0 * n0 * d1 * d1 / d0_3;
    } else {
      const double d0_2 = d0 * d0;
      const double d0_3 = d0_2 * d0;
      const double d0_4 = d0_3 * d0;
      out[i] = n3 / d0 -
        3.0 * n2 * d1 / d0_2 +
        3.0 * n1 * (2.0 * d1 * d1 / d0_3 - d2 / d0_2) +
        n0 * (-d3 / d0_2 + 6.0 * d1 * d2 / d0_3 -
              6.0 * d1 * d1 * d1 / d0_4);
    }
  }

  return out;
}

inline NumericMatrix gauss_kde_all_eval(const NumericVector& x,
                                        const NumericVector& x_obs,
                                        const NumericVector& y_obs,
                                        const double bw,
                                        const bool use_truncation = true,
                                        const double cutoff = DEFAULT_KDE_CUTOFF) {
  validate_kde_inputs(x, x_obs, y_obs, bw);

  const bool sorted_obs = is_strictly_increasing(x_obs);
  const bool truncate = use_truncation && sorted_obs;
  const R_xlen_t n = x.size();
  const R_xlen_t n_obs = x_obs.size();
  NumericMatrix out(n, 4);
  colnames(out) = CharacterVector::create("kde", "d1", "d2", "d3");
  const double denom_floor = std::numeric_limits<double>::min();

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];

    R_xlen_t j_lo = 0;
    R_xlen_t j_hi = n_obs - 1;

    if (truncate) {
      const std::pair<R_xlen_t, R_xlen_t> bounds =
        kde_support_bounds(xi, x_obs, cutoff * bw);
      j_lo = bounds.first;
      j_hi = bounds.second;
      if (j_lo > j_hi) {
        out(i, 0) = NA_REAL;
        out(i, 1) = NA_REAL;
        out(i, 2) = NA_REAL;
        out(i, 3) = NA_REAL;
        continue;
      }
    }

    double n0 = 0.0, d0 = 0.0;
    double n1 = 0.0, d1 = 0.0;
    double n2 = 0.0, d2 = 0.0;
    double n3 = 0.0, d3 = 0.0;

    kde_accumulate(xi, x_obs, y_obs, bw, j_lo, j_hi,
                   n0, d0, n1, d1, n2, d2, n3, d3);

    if (d0 <= denom_floor) {
      out(i, 0) = NA_REAL;
      out(i, 1) = NA_REAL;
      out(i, 2) = NA_REAL;
      out(i, 3) = NA_REAL;
      continue;
    }

    const double d0_2 = d0 * d0;
    const double d0_3 = d0_2 * d0;
    const double d0_4 = d0_3 * d0;

    out(i, 0) = n0 / d0;
    out(i, 1) = (n1 * d0 - n0 * d1) / d0_2;
    out(i, 2) = n2 / d0 - n0 * d2 / d0_2 - 2.0 * n1 * d1 / d0_2 +
      2.0 * n0 * d1 * d1 / d0_3;
    out(i, 3) = n3 / d0 -
      3.0 * n2 * d1 / d0_2 +
      3.0 * n1 * (2.0 * d1 * d1 / d0_3 - d2 / d0_2) +
      n0 * (-d3 / d0_2 + 6.0 * d1 * d2 / d0_3 -
            6.0 * d1 * d1 * d1 / d0_4);
  }

  return out;
}

inline NumericVector find_roots_on_grid(const NumericVector& x,
                                        const NumericVector& y_grid,
                                        const NumericVector& x_obs,
                                        const NumericVector& y_obs,
                                        const double bw,
                                        const double yi,
                                        const double tol,
                                        const int max_iter,
                                        const bool use_truncation,
                                        const double cutoff) {
  validate_root_inputs(x, yi, tol, max_iter);
  validate_numeric_vector(y_grid, "y_grid");

  if (y_grid.size() != x.size()) {
    stop("y_grid must have the same length as x.");
  }

  const R_xlen_t n = x.size();
  std::vector<double> roots;
  roots.reserve(static_cast<std::size_t>(n / 4 + 1));

  for (R_xlen_t i = 0; i < n; ++i) {
    const double fi = y_grid[i] - yi;
    if (std::abs(fi) <= tol &&
        !close_to_existing_root(roots, x[i], tol)) {
      roots.push_back(x[i]);
    }
  }

  for (R_xlen_t i = 0; i < n - 1; ++i) {
    const double f_left = y_grid[i] - yi;
    const double f_right = y_grid[i + 1] - yi;

    if (f_left == 0.0 || f_right == 0.0) continue;
    if (f_left * f_right < 0.0) {
      const double root = bisection_root_kde_deriv(
        x_obs, y_obs, bw, 1, use_truncation, cutoff,
        x[i], x[i + 1], yi, tol, max_iter
      );
      if (!close_to_existing_root(roots, root, tol)) {
        roots.push_back(root);
      }
    }
  }

  std::sort(roots.begin(), roots.end());
  return wrap(roots);
}

inline List find_extrema_from_grid(const NumericVector& x,
                                   const NumericVector& d1_grid,
                                   const NumericVector& x_obs,
                                   const NumericVector& y_obs,
                                   const double bw,
                                   const double tol,
                                   const int max_iter,
                                   const bool use_truncation,
                                   const double cutoff) {
  NumericVector x_roots = find_roots_on_grid(
    x, d1_grid, x_obs, y_obs, bw, 0.0, tol, max_iter, use_truncation, cutoff
  );

  if (x_roots.size() == 0) {
    return List::create(
      _["local_min"] = NumericVector::create(NA_REAL),
      _["local_max"] = NumericVector::create(NA_REAL),
      _["plateau"] = NumericVector::create(NA_REAL)
    );
  }

  const R_xlen_t n_roots = x_roots.size();
  NumericVector d2_roots(n_roots);

  for (R_xlen_t i = 0; i < n_roots; ++i) {
    d2_roots[i] = kde_deriv_scalar(
      x_roots[i], x_obs, y_obs, bw, 2, use_truncation, cutoff
    );
    if (!is_finite_double(d2_roots[i])) {
      stop("Second KDE derivative returned non-finite values at root positions.");
    }
  }

  std::vector<double> local_max;
  std::vector<double> local_min;
  std::vector<double> plateau;

  local_max.reserve(n_roots);
  local_min.reserve(n_roots);
  plateau.reserve(n_roots);

  for (R_xlen_t i = 0; i < n_roots; ++i) {
    const double d2_i = d2_roots[i];
    if (d2_i < -tol) {
      local_max.push_back(x_roots[i]);
    } else if (d2_i > tol) {
      local_min.push_back(x_roots[i]);
    } else {
      plateau.push_back(x_roots[i]);
    }
  }

  NumericVector out_local_max;
  if (local_max.empty()) {
    out_local_max = NumericVector::create(NA_REAL);
  } else {
    out_local_max = wrap(local_max);
  }

  NumericVector out_local_min;
  if (local_min.empty()) {
    out_local_min = NumericVector::create(NA_REAL);
  } else {
    out_local_min = wrap(local_min);
  }

  NumericVector out_plateau;
  if (plateau.empty()) {
    out_plateau = NumericVector::create(NA_REAL);
  } else {
    out_plateau = wrap(plateau);
  }

  return List::create(
    _["local_min"] = out_local_min,
    _["local_max"] = out_local_max,
    _["plateau"] = out_plateau
  );
}

inline NumericVector curvature_from_derivs(const NumericVector& d1,
                                           const NumericVector& d2,
                                           const bool absolute) {
  validate_numeric_vector(d1, "d1");
  validate_numeric_vector(d2, "d2");

  if (d1.size() != d2.size()) {
    stop("d1 and d2 must have the same length.");
  }

  const R_xlen_t n = d1.size();
  NumericVector out(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double d1_i = d1[i];
    const double d2_i = d2[i];
    const double denom = std::pow(1.0 + d1_i * d1_i, 1.5);
    const double numer = absolute ? std::abs(d2_i) : d2_i;
    out[i] = numer / denom;
  }

  return out;
}

} // namespace spectrum_math

// [[Rcpp::export]]
NumericVector cpp_gaussKDE(NumericVector x,
                           NumericVector x_obs,
                           NumericVector y_obs,
                           double bw) {
  return spectrum_math::gauss_kde_eval(x, x_obs, y_obs, bw, 0);
}

// [[Rcpp::export]]
NumericVector cpp_gaussKDE_1st_deriv(NumericVector x,
                                     NumericVector x_obs,
                                     NumericVector y_obs,
                                     double bw) {
  return spectrum_math::gauss_kde_eval(x, x_obs, y_obs, bw, 1);
}

// [[Rcpp::export]]
NumericVector cpp_gaussKDE_2nd_deriv(NumericVector x,
                                     NumericVector x_obs,
                                     NumericVector y_obs,
                                     double bw) {
  return spectrum_math::gauss_kde_eval(x, x_obs, y_obs, bw, 2);
}

// [[Rcpp::export]]
NumericVector cpp_gaussKDE_3rd_deriv(NumericVector x,
                                     NumericVector x_obs,
                                     NumericVector y_obs,
                                     double bw) {
  return spectrum_math::gauss_kde_eval(x, x_obs, y_obs, bw, 3);
}

// [[Rcpp::export]]
NumericMatrix cpp_gaussKDE_all(NumericVector x,
                               NumericVector x_obs,
                               NumericVector y_obs,
                               double bw) {
  return spectrum_math::gauss_kde_all_eval(x, x_obs, y_obs, bw);
}

// [[Rcpp::export]]
double cpp_bisection_root(Function fun,
                          double x_left,
                          double x_right,
                          double yi,
                          double tol = 1e-5,
                          int max_iter = 100) {
  if (!R_finite(x_left) || !R_finite(x_right)) {
    stop("x_left and x_right must be finite numbers.");
  }
  if (x_left == x_right) {
    stop("x_left and x_right must be different.");
  }
  if (!R_finite(yi)) {
    stop("yi must be a finite number.");
  }
  if (!R_finite(tol) || tol <= 0.0) {
    stop("tol must be a finite positive number.");
  }
  if (max_iter <= 0) {
    stop("max_iter must be a positive integer.");
  }

  return spectrum_math::bisection_root_cached(fun, x_left, x_right, yi, tol, max_iter);
}

// [[Rcpp::export]]
NumericVector cpp_find_roots(Function fun,
                             NumericVector x,
                             double yi,
                             double tol = 1e-5,
                             int max_iter = 100) {
  spectrum_math::validate_root_inputs(x, yi, tol, max_iter);

  NumericVector y = spectrum_math::eval_vector_function(fun, x);
  const R_xlen_t n = x.size();
  std::vector<double> roots;
  roots.reserve(static_cast<std::size_t>(n / 4 + 1));

  for (R_xlen_t i = 0; i < n; ++i) {
    const double fi = y[i] - yi;
    if (std::abs(fi) <= tol &&
        !spectrum_math::close_to_existing_root(roots, x[i], tol)) {
      roots.push_back(x[i]);
    }
  }

  for (R_xlen_t i = 0; i < n - 1; ++i) {
    const double f_left = y[i] - yi;
    const double f_right = y[i + 1] - yi;

    if (f_left == 0.0 || f_right == 0.0) continue;
    if (f_left * f_right < 0.0) {
      const double root = spectrum_math::bisection_root_cached(
        fun, x[i], x[i + 1], yi, tol, max_iter
      );
      if (!spectrum_math::close_to_existing_root(roots, root, tol)) {
        roots.push_back(root);
      }
    }
  }

  std::sort(roots.begin(), roots.end());
  return wrap(roots);
}

// [[Rcpp::export]]
NumericVector cpp_find_roots_on_grid(NumericVector x,
                                     NumericVector y_grid,
                                     NumericVector x_obs,
                                     NumericVector y_obs,
                                     double bw,
                                     double yi = 0.0,
                                     double tol = 1e-5,
                                     int max_iter = 100) {
  const bool use_truncation = spectrum_math::is_strictly_increasing(x_obs);
  return spectrum_math::find_roots_on_grid(
    x, y_grid, x_obs, y_obs, bw, yi, tol, max_iter,
    use_truncation, spectrum_math::DEFAULT_KDE_CUTOFF
  );
}

// [[Rcpp::export]]
List cpp_find_extrema_from_grid(NumericVector x,
                                NumericVector d1_grid,
                                NumericVector x_obs,
                                NumericVector y_obs,
                                double bw,
                                double tol = 1e-5,
                                int max_iter = 100) {
  const bool use_truncation = spectrum_math::is_strictly_increasing(x_obs);
  return spectrum_math::find_extrema_from_grid(
    x, d1_grid, x_obs, y_obs, bw, tol, max_iter,
    use_truncation, spectrum_math::DEFAULT_KDE_CUTOFF
  );
}

// [[Rcpp::export]]
NumericVector cpp_curvature(NumericVector x,
                            Function first_deriv,
                            Function second_deriv,
                            bool absolute = true) {
  spectrum_math::validate_numeric_vector(x, "x", true);

  NumericVector d1 = spectrum_math::eval_vector_function(first_deriv, x);
  NumericVector d2 = spectrum_math::eval_vector_function(second_deriv, x);

  return spectrum_math::curvature_from_derivs(d1, d2, absolute);
}

// [[Rcpp::export]]
NumericVector cpp_curvature_from_derivs(NumericVector d1,
                                      NumericVector d2,
                                      bool absolute = true) {
  return spectrum_math::curvature_from_derivs(d1, d2, absolute);
}
