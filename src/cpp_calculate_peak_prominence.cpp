#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_calculate_peak_prominence(NumericMatrix data,
                                            NumericMatrix peaks,
                                            std::string valley_type,
                                            double zero_tol) {
  
  // ----------------------------------------
  // Validate data
  // ----------------------------------------
  if (data.ncol() < 2) {
    stop("'data' must have at least two columns.");
  }
  
  if (data.nrow() < 3) {
    stop("'data' must contain at least 3 rows.");
  }
  
  if (peaks.ncol() < 2) {
    stop("'peaks' must have at least two columns.");
  }
  
  const int n = data.nrow();
  const int n_peaks = peaks.nrow();
  
  if (n_peaks == 0) {
    return NumericVector(0);
  }
  
  if (zero_tol < 0.0 || !R_finite(zero_tol)) {
    stop("'zero_tol' must be a single non-negative finite numeric value.");
  }
  
  // ----------------------------------------
  // Validate valley_type
  // ----------------------------------------
  if (valley_type != "higher" && valley_type != "lower") {
    stop("'valley_type' must be either 'higher' or 'lower'.");
  }
  
  // ----------------------------------------
  // Extract x and y
  // ----------------------------------------
  NumericVector x(n);
  NumericVector y(n);
  
  for (int i = 0; i < n; ++i) {
    x[i] = data(i, 0);
    y[i] = data(i, 1);
    
    if (!R_finite(x[i]) || !R_finite(y[i])) {
      stop("The first two columns of 'data' must contain finite numeric values.");
    }
  }
  
  // ----------------------------------------
  // Extract peak x and peak y
  // ----------------------------------------
  NumericVector x_peaks(n_peaks);
  NumericVector y_peaks(n_peaks);
  
  for (int i = 0; i < n_peaks; ++i) {
    x_peaks[i] = peaks(i, 0);
    y_peaks[i] = peaks(i, 1);
    
    if (!R_finite(x_peaks[i]) || !R_finite(y_peaks[i])) {
      stop("The first two columns of 'peaks' must contain finite numeric values.");
    }
  }
  
  // ----------------------------------------
  // Assumption:
  // data and peaks are already sorted by x in the R wrapper.
  // ----------------------------------------
  const double x_min = x[0];
  const double x_max = x[n - 1];
  
  NumericVector prominence(n_peaks);
  
  // ----------------------------------------
  // Calculate prominence for each peak
  // ----------------------------------------
  for (int i = 0; i < n_peaks; ++i) {
    
    const double xi_peak = x_peaks[i];
    const double yi_peak = y_peaks[i];
    
    // ----------------------------------------
    // Find left boundary:
    // nearest higher peak on the left
    // ----------------------------------------
    double x_left_bound = x_min;
    bool has_left_higher = false;
    
    for (int j = 0; j < n_peaks; ++j) {
      if (y_peaks[j] > yi_peak && x_peaks[j] < xi_peak) {
        if (!has_left_higher || x_peaks[j] > x_left_bound) {
          x_left_bound = x_peaks[j];
          has_left_higher = true;
        }
      }
    }
    
    // ----------------------------------------
    // Find right boundary:
    // nearest higher peak on the right
    // ----------------------------------------
    double x_right_bound = x_max;
    bool has_right_higher = false;
    
    for (int j = 0; j < n_peaks; ++j) {
      if (y_peaks[j] > yi_peak && x_peaks[j] > xi_peak) {
        if (!has_right_higher || x_peaks[j] < x_right_bound) {
          x_right_bound = x_peaks[j];
          has_right_higher = true;
        }
      }
    }
    
    // ----------------------------------------
    // Find left valley
    // ----------------------------------------
    double y_left_valley = yi_peak;
    bool has_left_point = false;
    
    for (int k = 0; k < n; ++k) {
      if (x[k] >= x_left_bound && x[k] <= xi_peak) {
        if (!has_left_point || y[k] < y_left_valley) {
          y_left_valley = y[k];
          has_left_point = true;
        }
      }
    }
    
    // ----------------------------------------
    // Find right valley
    // ----------------------------------------
    double y_right_valley = yi_peak;
    bool has_right_point = false;
    
    for (int k = 0; k < n; ++k) {
      if (x[k] >= xi_peak && x[k] <= x_right_bound) {
        if (!has_right_point || y[k] < y_right_valley) {
          y_right_valley = y[k];
          has_right_point = true;
        }
      }
    }
    
    // ----------------------------------------
    // Calculate prominence according to valley_type
    //
    // valley_type = "higher":
    //   Use the higher of the two valleys.
    //   This is conservative and corresponds to standard peak prominence.
    //
    // valley_type = "lower":
    //   Use the lower of the two valleys.
    //   This gives a larger prominence value.
    // ----------------------------------------
    double reference_valley;
    
    if (valley_type == "higher") {
      reference_valley = std::max(y_left_valley, y_right_valley);
    } else {
      reference_valley = std::min(y_left_valley, y_right_valley);
    }
    
    double prom = yi_peak - reference_valley;
    
    if (std::abs(prom) <= zero_tol) {
      prom = 0.0;
    }
    
    if (prom < 0.0) {
      prom = 0.0;
    }
    
    prominence[i] = prom;
  }
  
  return prominence;
}