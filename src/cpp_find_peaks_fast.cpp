#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// cpp_find_peaks_fast
// -----------------------------------------------------------------------------
// Finds local maxima from a two-column spectrum and keeps only the highest
// local maximum within each +/- hws_peaks window.
//
// plateau mode:
//   0 = middle
//   1 = first
//   2 = last
//
// The input data is assumed to be already validated in the R wrapper:
//   - two-column numeric matrix
//   - first column: x
//   - second column: y
//   - x values strictly increasing
//   - no NA/NaN/Inf
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix cpp_find_peaks_fast(NumericMatrix data,
                                  double hws_peaks = 10.0,
                                  int plateau = 0) {
  
  // ----------------------------------------
  // Basic validation
  // ----------------------------------------
  if (data.ncol() < 2) {
    stop("'data' must have at least two columns.");
  }
  
  const int n = data.nrow();
  
  if (n < 3) {
    stop("'data' must contain at least 3 rows.");
  }
  
  if (!R_finite(hws_peaks) || hws_peaks <= 0.0) {
    stop("'hws_peaks' must be a positive finite numeric value.");
  }
  
  if (plateau < 0 || plateau > 2) {
    stop("'plateau' must be 0, 1, or 2.");
  }
  
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
  // Find local maxima, including plateaus
  // ----------------------------------------
  std::vector<int> localmax_idx;
  localmax_idx.reserve(n / 10 + 1);
  
  int i = 1;
  
  while (i < n - 1) {
    
    // Strict local maximum
    if (y[i] > y[i - 1] && y[i] > y[i + 1]) {
      localmax_idx.push_back(i);
      ++i;
      continue;
    }
    
    // Plateau candidate: y[i] == y[i + 1]
    if (y[i] > y[i - 1] && y[i] == y[i + 1]) {
      
      const int start = i;
      int end = i + 1;
      
      while (end + 1 < n && y[end] == y[end + 1]) {
        ++end;
      }
      
      // A plateau is a local maximum only if it descends after the plateau.
      if (end < n - 1 && y[end] > y[end + 1]) {
        int idx_selected;
        
        if (plateau == 1) {
          // first
          idx_selected = start;
        } else if (plateau == 2) {
          // last
          idx_selected = end;
        } else {
          // middle
          idx_selected = start + (end - start) / 2;
        }
        
        localmax_idx.push_back(idx_selected);
      }
      
      i = end + 1;
      continue;
    }
    
    ++i;
  }
  
  const int n_localmax = static_cast<int>(localmax_idx.size());
  
  if (n_localmax == 0) {
    NumericMatrix out(0, 2);
    colnames(out) = CharacterVector::create("x", "y");
    return out;
  }
  
  NumericVector localmax_x(n_localmax);
  NumericVector localmax_y(n_localmax);
  
  for (int j = 0; j < n_localmax; ++j) {
    const int idx = localmax_idx[j];
    localmax_x[j] = x[idx];
    localmax_y[j] = y[idx];
  }
  
  // ----------------------------------------
  // Keep only the highest local maximum within each window
  // Same tie behavior as the R version:
  // if multiple local maxima have the same max y in a window, keep all of them.
  // ----------------------------------------
  std::vector<int> peak_idx;
  peak_idx.reserve(n_localmax);
  
  for (int a = 0; a < n_localmax; ++a) {
    
    const double left = localmax_x[a] - hws_peaks;
    const double right = localmax_x[a] + hws_peaks;
    
    double max_y = R_NegInf;
    
    for (int b = 0; b < n_localmax; ++b) {
      if (localmax_x[b] >= left && localmax_x[b] <= right) {
        if (localmax_y[b] > max_y) {
          max_y = localmax_y[b];
        }
      }
    }
    
    if (localmax_y[a] == max_y) {
      peak_idx.push_back(a);
    }
  }
  
  const int n_peaks = static_cast<int>(peak_idx.size());
  NumericMatrix out(n_peaks, 2);
  
  for (int p = 0; p < n_peaks; ++p) {
    const int idx = peak_idx[p];
    out(p, 0) = localmax_x[idx];
    out(p, 1) = localmax_y[idx];
  }
  
  colnames(out) = CharacterVector::create("x", "y");
  return out;
}
