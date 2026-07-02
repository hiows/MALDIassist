#include <Rcpp.h>
#include <map>
#include <utility>
using namespace Rcpp;

struct SavgolCoefKey {
  int hws;
  int pno;

  bool operator<(const SavgolCoefKey& other) const {
    if (hws != other.hws) return hws < other.hws;
    return pno < other.pno;
  }
};

static std::map<SavgolCoefKey, std::vector<std::vector<double> > > savgol_coef_cache;

// -----------------------------------------------------------------------------
// Solve linear system A * x = b using Gauss-Jordan elimination.
// This small helper is used internally to compute Savitzky-Golay coefficients.
// -----------------------------------------------------------------------------
std::vector<double> solve_linear_system(std::vector<std::vector<double> > A,
                                        std::vector<double> b) {
  const int n = A.size();

  for (int i = 0; i < n; ++i) {

    // Find pivot
    int pivot = i;
    double max_abs = std::abs(A[i][i]);

    for (int r = i + 1; r < n; ++r) {
      double val = std::abs(A[r][i]);
      if (val > max_abs) {
        max_abs = val;
        pivot = r;
      }
    }

    if (max_abs == 0.0) {
      stop("Singular matrix encountered while calculating Savitzky-Golay coefficients.");
    }

    // Swap rows if needed
    if (pivot != i) {
      std::swap(A[i], A[pivot]);
      std::swap(b[i], b[pivot]);
    }

    // Normalize pivot row
    double diag = A[i][i];
    for (int c = 0; c < n; ++c) {
      A[i][c] /= diag;
    }
    b[i] /= diag;

    // Eliminate other rows
    for (int r = 0; r < n; ++r) {
      if (r == i) continue;

      double factor = A[r][i];
      for (int c = 0; c < n; ++c) {
        A[r][c] -= factor * A[i][c];
      }
      b[r] -= factor * b[i];
    }
  }

  return b;
}

// -----------------------------------------------------------------------------
// Calculate one row of Savitzky-Golay smoothing coefficients.
// The returned row corresponds to the intercept estimate at the target position.
// -----------------------------------------------------------------------------
std::vector<double> calculate_savgol_coef_row(const int window_size,
                                              const int pno,
                                              const int target_index) {
  const int n_coef = pno + 1;

  std::vector<std::vector<double> > XtX(
    n_coef,
    std::vector<double>(n_coef, 0.0)
  );

  std::vector<std::vector<double> > X(
    window_size,
    std::vector<double>(n_coef, 0.0)
  );

  // Same offset logic as R code:
  // offset <- seq_len(window_size) - i
  // where target_index is zero-based and i is one-based in R.
  for (int r = 0; r < window_size; ++r) {
    double offset = static_cast<double>((r + 1) - (target_index + 1));
    double value = 1.0;

    for (int c = 0; c < n_coef; ++c) {
      X[r][c] = value;
      value *= offset;
    }
  }

  // XtX = t(X) %*% X
  for (int i = 0; i < n_coef; ++i) {
    for (int j = 0; j < n_coef; ++j) {
      double sum = 0.0;
      for (int r = 0; r < window_size; ++r) {
        sum += X[r][i] * X[r][j];
      }
      XtX[i][j] = sum;
    }
  }

  // Coefficient row = first row of solve(XtX, t(X)).
  // Equivalent to solving XtX * beta = e0, then coef = t(beta) %*% t(X).
  std::vector<double> e0(n_coef, 0.0);
  e0[0] = 1.0;

  std::vector<double> beta = solve_linear_system(XtX, e0);

  std::vector<double> coef(window_size, 0.0);
  for (int r = 0; r < window_size; ++r) {
    double sum = 0.0;
    for (int c = 0; c < n_coef; ++c) {
      sum += beta[c] * X[r][c];
    }
    coef[r] = sum;
  }

  return coef;
}

// [[Rcpp::export]]
NumericVector cpp_savitzkyGolay_filter(NumericVector y,
                                       int hws,
                                       int pno) {

  // ----------------------------------------
  // Input validation
  // ----------------------------------------
  const int n = y.size();

  if (n < 1) {
    stop("'y' must contain at least one value.");
  }

  if (hws < 1) {
    stop("'hws' must be greater than or equal to 1.");
  }

  if (pno < 0) {
    stop("'pno' must be greater than or equal to 0.");
  }

  const int window_size = 2 * hws + 1;

  if (window_size <= pno) {
    stop("The full window size must be larger than the polynomial order.");
  }

  if (n < window_size) {
    stop("The length of 'y' must be greater than or equal to the full window size.");
  }

  for (int i = 0; i < n; ++i) {
    if (!R_finite(y[i])) {
      stop("'y' must contain only finite non-missing values.");
    }
  }

  // ----------------------------------------
  // Calculate or retrieve cached Savitzky-Golay coefficients
  // ----------------------------------------
  const SavgolCoefKey cache_key = {hws, pno};
  std::vector<std::vector<double> > coef;

  std::map<SavgolCoefKey, std::vector<std::vector<double> > >::const_iterator cache_it =
    savgol_coef_cache.find(cache_key);

  if (cache_it != savgol_coef_cache.end()) {
    coef = cache_it->second;
  } else {
    coef.assign(
      window_size,
      std::vector<double>(window_size, 0.0)
    );

    for (int i = 0; i <= hws; ++i) {
      coef[i] = calculate_savgol_coef_row(window_size, pno, i);
    }

    for (int i = hws + 1; i < window_size; ++i) {
      int mirror_row = window_size - 1 - i;
      for (int j = 0; j < window_size; ++j) {
        coef[i][j] = coef[mirror_row][window_size - 1 - j];
      }
    }

    savgol_coef_cache[cache_key] = coef;
  }

  NumericVector out(n);

  // ----------------------------------------
  // Apply left boundary smoothing
  // ----------------------------------------
  for (int i = 0; i < hws; ++i) {
    double sum = 0.0;
    for (int j = 0; j < window_size; ++j) {
      sum += coef[i][j] * y[j];
    }
    out[i] = sum;
  }

  // ----------------------------------------
  // Apply center smoothing
  // ----------------------------------------
  const int center_row = hws;

  for (int i = hws; i <= n - hws - 1; ++i) {
    double sum = 0.0;
    int start = i - hws;

    for (int j = 0; j < window_size; ++j) {
      sum += coef[center_row][j] * y[start + j];
    }

    out[i] = sum;
  }

  // ----------------------------------------
  // Apply right boundary smoothing
  // ----------------------------------------
  for (int i = n - hws; i < n; ++i) {
    int coef_row = hws + 1 + (i - (n - hws));
    double sum = 0.0;
    int start = n - window_size;

    for (int j = 0; j < window_size; ++j) {
      sum += coef[coef_row][j] * y[start + j];
    }

    out[i] = sum;
  }

  return out;
}
