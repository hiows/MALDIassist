#include <Rcpp.h>
#include <deque>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ------------------------------------------------------------
// Utility: validate a numeric vector for baseline correction
// ------------------------------------------------------------
static void validate_intensity_vector(const NumericVector& x, const char* name) {
  if (x.size() == 0) {
    stop("'%s' must not be empty.", name);
  }

  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i])) {
      stop("'%s' contains NA, NaN, Inf, or -Inf values.", name);
    }
  }
}

// ------------------------------------------------------------
// Utility: moving minimum using a deque
//
// This computes a centered moving minimum with boundary clipping.
// For each index i, it uses the interval:
//   [max(0, i - half_window), min(n - 1, i + half_window)]
//
// Complexity: O(n)
// ------------------------------------------------------------
static NumericVector moving_min_centered(const NumericVector& x, const int half_window) {
  const int n = x.size();
  NumericVector out(n);

  std::deque<int> dq;
  int right_added = -1;

  for (int i = 0; i < n; ++i) {
    const int left = std::max(0, i - half_window);
    const int right = std::min(n - 1, i + half_window);

    // Add new right-side elements into deque.
    while (right_added < right) {
      ++right_added;

      while (!dq.empty() && x[dq.back()] >= x[right_added]) {
        dq.pop_back();
      }
      dq.push_back(right_added);
    }

    // Remove elements outside the left boundary.
    while (!dq.empty() && dq.front() < left) {
      dq.pop_front();
    }

    out[i] = x[dq.front()];
  }

  return out;
}

// ------------------------------------------------------------
// Utility: moving maximum using a deque
//
// This computes a centered moving maximum with boundary clipping.
// For each index i, it uses the interval:
//   [max(0, i - half_window), min(n - 1, i + half_window)]
//
// Complexity: O(n)
// ------------------------------------------------------------
static NumericVector moving_max_centered(const NumericVector& x, const int half_window) {
  const int n = x.size();
  NumericVector out(n);

  std::deque<int> dq;
  int right_added = -1;

  for (int i = 0; i < n; ++i) {
    const int left = std::max(0, i - half_window);
    const int right = std::min(n - 1, i + half_window);

    // Add new right-side elements into deque.
    while (right_added < right) {
      ++right_added;

      while (!dq.empty() && x[dq.back()] <= x[right_added]) {
        dq.pop_back();
      }
      dq.push_back(right_added);
    }

    // Remove elements outside the left boundary.
    while (!dq.empty() && dq.front() < left) {
      dq.pop_front();
    }

    out[i] = x[dq.front()];
  }

  return out;
}

// ------------------------------------------------------------
// cpp_SNIP
//
// Sensitive Nonlinear Iterative Peak-clipping baseline estimation.
//
// This function estimates the baseline from a one-dimensional intensity
// vector. The default implementation follows the common SNIP-style update:
//
//   b_i <- min(b_i, (b_{i-k} + b_{i+k}) / 2)
//
// for decreasing or increasing clipping window sizes. For MALDI-TOF spectra,
// decreasing = true is often used in practice after optional smoothing.
//
// Parameters
// ----------
// intensity : NumericVector
//   Spectrum intensity values.
// iterations : int
//   Number of clipping iterations. Common values are around 20-100.
// decreasing : bool
//   If true, window size decreases from iterations to 1.
//   If false, window size increases from 1 to iterations.
// nonnegative : bool
//   If true, negative baseline values are clipped to zero.
//
// Returns
// -------
// NumericVector
//   Estimated baseline with the same length as intensity.
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_SNIP(NumericVector intensity,
                       int iterations = 50,
                       bool decreasing = true,
                       bool nonnegative = true) {

  validate_intensity_vector(intensity, "intensity");

  if (iterations < 1) {
    stop("'iterations' must be a positive integer.");
  }

  const int n = intensity.size();
  NumericVector baseline = clone(intensity);
  NumericVector previous(n);

  // If the vector is too short, there is no meaningful clipping interval.
  if (n < 3) {
    if (nonnegative) {
      for (int i = 0; i < n; ++i) {
        if (baseline[i] < 0.0) baseline[i] = 0.0;
      }
    }
    return baseline;
  }

  // The maximum usable half-window is bounded by the data length.
  const int max_k = std::min(iterations, (n - 1) / 2);

  if (decreasing) {
    for (int k = max_k; k >= 1; --k) {
      std::swap(baseline, previous);

      for (int i = k; i < n - k; ++i) {
        const double candidate = 0.5 * (previous[i - k] + previous[i + k]);
        if (candidate < previous[i]) {
          baseline[i] = candidate;
        } else {
          baseline[i] = previous[i];
        }
      }

      for (int i = 0; i < k; ++i) {
        baseline[i] = previous[i];
      }
      for (int i = n - k; i < n; ++i) {
        baseline[i] = previous[i];
      }
    }
  } else {
    for (int k = 1; k <= max_k; ++k) {
      std::swap(baseline, previous);

      for (int i = k; i < n - k; ++i) {
        const double candidate = 0.5 * (previous[i - k] + previous[i + k]);
        if (candidate < previous[i]) {
          baseline[i] = candidate;
        } else {
          baseline[i] = previous[i];
        }
      }

      for (int i = 0; i < k; ++i) {
        baseline[i] = previous[i];
      }
      for (int i = n - k; i < n; ++i) {
        baseline[i] = previous[i];
      }
    }
  }

  if (nonnegative) {
    for (int i = 0; i < n; ++i) {
      if (baseline[i] < 0.0) baseline[i] = 0.0;
    }
  }

  return baseline;
}

// ------------------------------------------------------------
// cpp_SNIP_spectrum
//
// Convenience wrapper for a two-column spectrum matrix.
// spectrum[, 1] = m/z
// spectrum[, 2] = intensity
//
// Returns a data.frame with:
//   mz, intensity, baseline, corrected
// ------------------------------------------------------------

// [[Rcpp::export]]
DataFrame cpp_SNIP_spectrum(NumericMatrix spectrum,
                            int iterations = 50,
                            bool decreasing = true,
                            bool nonnegative = true,
                            bool corrected_nonnegative = true) {

  if (spectrum.ncol() < 2) {
    stop("'spectrum' must have at least two columns: m/z and intensity.");
  }

  const int n = spectrum.nrow();
  if (n == 0) {
    stop("'spectrum' must not be empty.");
  }

  NumericVector mz(n);
  NumericVector intensity(n);

  for (int i = 0; i < n; ++i) {
    mz[i] = spectrum(i, 0);
    intensity[i] = spectrum(i, 1);

    if (!R_finite(mz[i]) || !R_finite(intensity[i])) {
      stop("'spectrum' contains NA, NaN, Inf, or -Inf values.");
    }
  }

  NumericVector baseline = cpp_SNIP(intensity, iterations, decreasing, nonnegative);
  NumericVector corrected(n);

  for (int i = 0; i < n; ++i) {
    corrected[i] = intensity[i] - baseline[i];
    if (corrected_nonnegative && corrected[i] < 0.0) {
      corrected[i] = 0.0;
    }
  }

  return DataFrame::create(
    Named("mz") = mz,
    Named("intensity") = intensity,
    Named("baseline") = baseline,
    Named("corrected") = corrected
  );
}

// ------------------------------------------------------------
// cpp_TopHat
//
// Morphological Top-Hat baseline estimation for 1D spectra.
//
// This function estimates a slowly varying baseline by morphological opening:
//
//   opening(x) = dilation(erosion(x))
//
// where erosion is a moving minimum and dilation is a moving maximum.
// The estimated baseline is the opened signal. Baseline-corrected intensity is:
//
//   corrected = intensity - baseline
//
// Here, half_window defines the structuring element radius.
// Full window size = 2 * half_window + 1.
//
// Parameters
// ----------
// intensity : NumericVector
//   Spectrum intensity values.
// half_window : int
//   Radius of the structuring element.
// nonnegative : bool
//   If true, negative baseline values are clipped to zero.
//
// Returns
// -------
// NumericVector
//   Estimated baseline with the same length as intensity.
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_TopHat(NumericVector intensity,
                         int half_window = 50,
                         bool nonnegative = true) {

  validate_intensity_vector(intensity, "intensity");

  if (half_window < 1) {
    stop("'half_window' must be a positive integer.");
  }

  const int n = intensity.size();
  const int hw = std::min(half_window, std::max(1, n - 1));

  NumericVector eroded = moving_min_centered(intensity, hw);
  NumericVector baseline = moving_max_centered(eroded, hw);

  if (nonnegative) {
    for (int i = 0; i < n; ++i) {
      if (baseline[i] < 0.0) baseline[i] = 0.0;
    }
  }

  return baseline;
}

// ------------------------------------------------------------
// cpp_TopHat_spectrum
//
// Convenience wrapper for a two-column spectrum matrix.
// spectrum[, 1] = m/z
// spectrum[, 2] = intensity
//
// Returns a data.frame with:
//   mz, intensity, baseline, corrected
// ------------------------------------------------------------

// [[Rcpp::export]]
DataFrame cpp_TopHat_spectrum(NumericMatrix spectrum,
                              int half_window = 50,
                              bool nonnegative = true,
                              bool corrected_nonnegative = true) {

  if (spectrum.ncol() < 2) {
    stop("'spectrum' must have at least two columns: m/z and intensity.");
  }

  const int n = spectrum.nrow();
  if (n == 0) {
    stop("'spectrum' must not be empty.");
  }

  NumericVector mz(n);
  NumericVector intensity(n);

  for (int i = 0; i < n; ++i) {
    mz[i] = spectrum(i, 0);
    intensity[i] = spectrum(i, 1);

    if (!R_finite(mz[i]) || !R_finite(intensity[i])) {
      stop("'spectrum' contains NA, NaN, Inf, or -Inf values.");
    }
  }

  NumericVector baseline = cpp_TopHat(intensity, half_window, nonnegative);
  NumericVector corrected(n);

  for (int i = 0; i < n; ++i) {
    corrected[i] = intensity[i] - baseline[i];
    if (corrected_nonnegative && corrected[i] < 0.0) {
      corrected[i] = 0.0;
    }
  }

  return DataFrame::create(
    Named("mz") = mz,
    Named("intensity") = intensity,
    Named("baseline") = baseline,
    Named("corrected") = corrected
  );
}
