#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double EPS = std::numeric_limits<double>::epsilon();

inline bool is_finite_double(double x) {
  return R_finite(x);
}

NumericVector extract_numeric_column(SEXP obj, int col_index_zero_based,
                                     const std::string& obj_name) {
  if (Rf_isMatrix(obj)) {
    NumericMatrix mat(obj);

    if (mat.ncol() <= col_index_zero_based) {
      stop("'%s' must have at least %d columns.",
           obj_name, col_index_zero_based + 1);
    }

    NumericVector out(mat.nrow());
    for (int i = 0; i < mat.nrow(); ++i) {
      out[i] = mat(i, col_index_zero_based);
    }
    return out;
  }

  if (Rf_inherits(obj, "data.frame")) {
    DataFrame df(obj);

    if (df.size() <= col_index_zero_based) {
      stop("'%s' must have at least %d columns.",
           obj_name, col_index_zero_based + 1);
    }

    Function as_numeric("as.numeric");
    return as<NumericVector>(as_numeric(df[col_index_zero_based]));
  }

  stop("'%s' must be a data frame or matrix.", obj_name);
  return NumericVector();
}

int object_nrow(SEXP obj, const std::string& obj_name) {
  if (Rf_isMatrix(obj)) {
    NumericMatrix mat(obj);
    return mat.nrow();
  }

  if (Rf_inherits(obj, "data.frame")) {
    DataFrame df(obj);
    if (df.size() == 0) return 0;
    return Rf_length(df[0]);
  }

  stop("'%s' must be a data frame or matrix.", obj_name);
  return 0;
}

void validate_no_na_finite(const NumericVector& x, const std::string& name) {
  for (int i = 0; i < x.size(); ++i) {
    if (NumericVector::is_na(x[i])) {
      stop("'%s' contains NA values.", name);
    }
    if (!is_finite_double(x[i])) {
      stop("'%s' contains non-finite values.", name);
    }
  }
}

bool has_duplicate_sorted_input(const NumericVector& x) {
  std::vector<double> tmp = as<std::vector<double> >(x);
  std::sort(tmp.begin(), tmp.end());

  for (std::size_t i = 1; i < tmp.size(); ++i) {
    if (tmp[i] == tmp[i - 1]) return true;
  }
  return false;
}

bool is_strictly_increasing(const NumericVector& x) {
  for (int i = 1; i < x.size(); ++i) {
    if (!(x[i] > x[i - 1])) return false;
  }
  return true;
}

std::vector<double> positive_values(const std::vector<double>& x) {
  std::vector<double> out;
  out.reserve(x.size());

  for (double v : x) {
    if (is_finite_double(v) && v > 0.0) out.push_back(v);
  }
  return out;
}

double median_cpp(std::vector<double> x) {
  if (x.empty()) return NA_REAL;

  std::sort(x.begin(), x.end());
  const std::size_t n = x.size();
  const std::size_t mid = n / 2;

  if (n % 2 == 1) {
    return x[mid];
  }

  return (x[mid - 1] + x[mid]) / 2.0;
}

double mad_cpp(const std::vector<double>& x, double constant = 1.4826) {
  if (x.empty()) return NA_REAL;

  const double med = median_cpp(x);
  if (!is_finite_double(med)) return NA_REAL;

  std::vector<double> abs_dev;
  abs_dev.reserve(x.size());

  for (double v : x) {
    abs_dev.push_back(std::abs(v - med));
  }

  const double mad_raw = median_cpp(abs_dev);
  if (!is_finite_double(mad_raw)) return NA_REAL;

  return constant * mad_raw;
}

// R quantile type 7, equivalent to stats::quantile(..., type = 7)
double quantile_type7_cpp(std::vector<double> x, double prob) {
  if (x.empty()) return NA_REAL;
  if (prob < 0.0 || prob > 1.0) return NA_REAL;

  std::sort(x.begin(), x.end());
  const std::size_t n = x.size();

  if (n == 1) return x[0];
  if (prob == 0.0) return x.front();
  if (prob == 1.0) return x.back();

  const double h = 1.0 + (static_cast<double>(n) - 1.0) * prob;
  const std::size_t j = static_cast<std::size_t>(std::floor(h));
  const double gamma = h - static_cast<double>(j);

  const double xj = x[j - 1];
  const double xj1 = x[std::min(j, n - 1)];

  return (1.0 - gamma) * xj + gamma * xj1;
}

double mean_cpp(const std::vector<double>& x) {
  if (x.empty()) return NA_REAL;

  double sum = 0.0;
  int n = 0;

  for (double v : x) {
    if (is_finite_double(v)) {
      sum += v;
      ++n;
    }
  }

  if (n == 0) return NA_REAL;
  return sum / static_cast<double>(n);
}

inline double transform_intensity(double y, const std::string& normalization) {
  if (normalization == "raw") {
    return y;
  }

  if (normalization == "sqrt") {
    return std::sqrt(y);
  }

  if (normalization == "log10") {
    return std::log10(1.0 + y);
  }

  stop("'normalization' must be one of 'raw', 'sqrt', or 'log10'.");
  return NA_REAL;
}

} // namespace

// cpp_estimate_peak_strength (C++ backend)
//
// Estimates a bounded peak strength score for each supplied peak. The signal
// intensity can be used as-is or transformed ("raw", "sqrt", "log10") before
// calculation. A robust global strength scale is estimated as
// median(y_global) + k * MAD(y_global); a local scale is estimated between the
// left and right valleys of each peak. The score is
// local / (global + local), bounded in [0, 1] when both scales are positive.
// Peaks for which local valleys cannot be defined return NA_real_.
//
// Internal roxygen documentation lives in R/auxiliary_cpp.R.
// [[Rcpp::export]]
NumericVector cpp_estimate_peak_strength(SEXP data,
                                         SEXP peaks,
                                         double k = 1.0,
                                         std::string normalization = "raw") {
  // ----------------------------------------
  // 1. Validate data
  // ----------------------------------------
  const int n_data = object_nrow(data, "data");

  if (n_data < 3) {
    stop("'data' must contain at least 3 rows.");
  }

  NumericVector x_raw = extract_numeric_column(data, 0, "data");
  NumericVector y_raw = extract_numeric_column(data, 1, "data");

  if (x_raw.size() != y_raw.size()) {
    stop("The first and second columns of 'data' must have the same length.");
  }

  validate_no_na_finite(x_raw, "x values in 'data'");
  validate_no_na_finite(y_raw, "y values in 'data'");

  if (has_duplicate_sorted_input(x_raw)) {
    stop("'x' values in 'data' must not be duplicated.");
  }

  if (!is_strictly_increasing(x_raw)) {
    stop("'x' values in 'data' must be strictly increasing.");
  }

  // ----------------------------------------
  // 2. Validate parameters
  // ----------------------------------------
  if (!is_finite_double(k) || NumericVector::is_na(k)) {
    stop("'k' must be a single finite numeric value.");
  }

  if (k < 0.0) {
    stop("'k' must be greater than or equal to 0.");
  }

  if (normalization != "raw" &&
      normalization != "sqrt" &&
      normalization != "log10") {
    stop("'normalization' must be one of 'raw', 'sqrt', or 'log10'.");
  }

  for (int i = 0; i < y_raw.size(); ++i) {
    if (y_raw[i] < 0.0) {
      stop("'y' values must be non-negative when normalization is 'raw', 'sqrt', or 'log10'.");
    }
  }

  // ----------------------------------------
  // 3. Transform intensity
  // ----------------------------------------
  const int n = x_raw.size();
  std::vector<double> x(n);
  std::vector<double> y(n);

  for (int i = 0; i < n; ++i) {
    x[i] = x_raw[i];
    y[i] = transform_intensity(y_raw[i], normalization);
  }

  // ----------------------------------------
  // 4. Estimate global strength scale
  // ----------------------------------------
  std::vector<double> yy_global = positive_values(y);

  if (yy_global.empty()) {
    return NumericVector(0);
  }

  double global_strength = median_cpp(yy_global) + mad_cpp(yy_global) * k;

  if (!is_finite_double(global_strength) || global_strength <= EPS) {
    global_strength = quantile_type7_cpp(yy_global, 0.95);
  }

  if (!is_finite_double(global_strength) || global_strength <= EPS) {
    global_strength = EPS;
  }

  // ----------------------------------------
  // 5. Validate and sort peaks
  // ----------------------------------------
  const int n_peaks = object_nrow(peaks, "peaks");

  if (n_peaks == 0) {
    return NumericVector(0);
  }

  NumericVector peak_x_vec = extract_numeric_column(peaks, 0, "peaks");
  NumericVector peak_y_vec = extract_numeric_column(peaks, 1, "peaks");

  if (peak_x_vec.size() != peak_y_vec.size()) {
    stop("The first and second columns of 'peaks' must have the same length.");
  }

  validate_no_na_finite(peak_x_vec, "peak x values in 'peaks'");
  validate_no_na_finite(peak_y_vec, "peak y values in 'peaks'");

  std::vector<std::pair<double, double> > peak_pairs;
  peak_pairs.reserve(peak_x_vec.size());

  for (int i = 0; i < peak_x_vec.size(); ++i) {
    peak_pairs.emplace_back(peak_x_vec[i], peak_y_vec[i]);
  }

  std::sort(
    peak_pairs.begin(),
    peak_pairs.end(),
    [](const std::pair<double, double>& a,
       const std::pair<double, double>& b) {
      return a.first < b.first;
    }
  );

  // Boundary points are used only for valley search.
  std::vector<double> boundary_x;
  boundary_x.reserve(peak_pairs.size() + 2);

  boundary_x.push_back(x.front());
  for (const auto& p : peak_pairs) {
    boundary_x.push_back(p.first);
  }
  boundary_x.push_back(x.back());

  const int n_with_boundary = boundary_x.size();
  NumericVector peak_strength(peak_pairs.size());

  // ----------------------------------------
  // 6. Estimate peak strength for each peak
  // ----------------------------------------
  for (int i = 1; i < n_with_boundary - 1; ++i) {
    const double xi = boundary_x[i];

    // ------------------------------------
    // 6-1. Find left valley
    // ------------------------------------
    bool has_left = false;
    double x_left_valley = NA_REAL;
    double min_left_y = std::numeric_limits<double>::infinity();

    for (int j = 0; j < n; ++j) {
      if (x[j] > boundary_x[i - 1] && x[j] < xi) {
        if (y[j] < min_left_y) {
          min_left_y = y[j];
          x_left_valley = x[j];
          has_left = true;
        }
      }
    }

    if (!has_left) {
      peak_strength[i - 1] = NA_REAL;
      continue;
    }

    // ------------------------------------
    // 6-2. Find right valley
    // ------------------------------------
    bool has_right = false;
    double x_right_valley = NA_REAL;
    double min_right_y = std::numeric_limits<double>::infinity();

    for (int j = 0; j < n; ++j) {
      if (x[j] > xi && x[j] < boundary_x[i + 1]) {
        if (y[j] < min_right_y) {
          min_right_y = y[j];
          x_right_valley = x[j];
          has_right = true;
        }
      }
    }

    if (!has_right) {
      peak_strength[i - 1] = NA_REAL;
      continue;
    }

    // ------------------------------------
    // 6-3. Estimate local strength scale
    // ------------------------------------
    std::vector<double> yy_local;

    for (int j = 0; j < n; ++j) {
      if (x[j] > x_left_valley &&
          x[j] < x_right_valley &&
          is_finite_double(y[j])) {
        yy_local.push_back(y[j]);
      }
    }

    if (yy_local.empty()) {
      peak_strength[i - 1] = NA_REAL;
      continue;
    }

    double local_strength = mean_cpp(yy_local) + mad_cpp(yy_local) * k;

    if (!is_finite_double(local_strength) || local_strength <= EPS) {
      std::vector<double> yy_local_pos = positive_values(yy_local);

      if (yy_local_pos.empty()) {
        local_strength = EPS;
      } else {
        local_strength = quantile_type7_cpp(yy_local_pos, 0.95);
      }
    }

    if (!is_finite_double(local_strength) || local_strength <= EPS) {
      local_strength = EPS;
    }

    // ------------------------------------
    // 6-4. Peak strength score
    // ------------------------------------
    const double strength_i =
      local_strength / (global_strength + local_strength);

    if (!is_finite_double(strength_i)) {
      peak_strength[i - 1] = NA_REAL;
    } else {
      peak_strength[i - 1] = strength_i;
    }
  }

  return peak_strength;
}