# MALDIassist

<!-- badges: start -->
<!-- badges: end -->

**MALDIassist** provides a set of mathematical utilities for **MALDI-TOF mass spectrometry** workflows in R. It covers the full path from raw Bruker spectra to a cohort-level peak matrix:

- Loading Bruker MALDI-TOF spectra
- Smoothing and baseline correction
- Gaussian KDE-based peak detection (including shoulder peaks)
- Peak-quality metrics and filtering (intensity / prominence / strength)
- Cohort feature analysis (frequent m/z discovery, matched peak matrix)

Performance-critical routines are implemented in C++ (via [Rcpp](https://www.rcpp.org/)), and list-based functions support multi-core parallel processing.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hiows/MALDIassist")
```

---

## Quick start

The example below mirrors the full pipeline: load → preprocess → detect peaks → filter → visualize.

### 0. Load the library

```r
library(MALDIassist)
```

### 1. Load spectra

Point `load_maldi_spectra()` at a directory of Bruker flex data files. It returns a named list of raw spectra.

```r
raw_spectra <- load_maldi_spectra(spectra_dir = "data/")
```

### 2. Preprocessing

Apply Savitzky-Golay smoothing and baseline subtraction. List-based functions accept `n_cores` for parallel processing.

```r
preprocessed_spectra <- preprocess_maldi_spectra(
  spectra       = raw_spectra,
  hws_sg        = 10,      # half-window size for Savitzky-Golay
  pno_sg        = 3,       # polynomial order
  baseline_type = "snip",  # baseline algorithm
  iter_snip     = 100,     # SNIP iterations
  n_cores       = 4
)
```

Compare a raw spectrum (gray) against its preprocessed version (red) over a chosen m/z window:

```r
example_range <- c(12000, 15000)
i    <- 1
spec <- raw_spectra[[i]]
pp   <- preprocessed_spectra[[i]]

plot(
  subset(spec, mz >= example_range[1] & mz <= example_range[2]),
  type = "l", col = "gray60", lwd = 1.5,
  ylim = c(0, max(spec[, 2])),
  xlab = "m/z", ylab = "Intensity", frame = FALSE
)
lines(
  subset(pp, mz >= example_range[1] & mz <= example_range[2]),
  col = "red", lwd = 2
)
```

### 3. Find and filter peaks

`find_peaks_spectra()` detects peaks (including shoulder peaks) using a Gaussian KDE approach; `filter_peaks_spectra()` removes low-quality peaks by intensity, prominence, and strength cutoffs.

```r
peaks_list <- find_peaks_spectra(
  spectra                    = preprocessed_spectra,
  bw                         = 1,      # KDE bandwidth
  hws_peaks                  = 10,
  weight_type                = "raw",
  cutoff_kappa_peak_strength = 0.3,
  n_cores                    = 4
)

filtered_peaks_list <- filter_peaks_spectra(
  spectra                 = preprocessed_spectra,
  peaks_list              = peaks_list,
  cutoff_peak_intensity   = 100,
  cutoff_peak_prominence  = 100,
  cutoff_peak_strength    = 0.5,
  normalization_type      = "raw"
)
```

Overlay the filtered peaks on the spectrum plot from step 2:

```r
filtered_peaks <- filtered_peaks_list[[i]]
segments(
  filtered_peaks$mz, 0,
  filtered_peaks$mz, filtered_peaks$intensity,
  col = "blue", lwd = 1.5
)
```

### 4. Visualization

```r
plot_spectra(spectra = preprocessed_spectra)
```

---

## Main functions

| Function | Purpose |
| --- | --- |
| `load_maldi_spectra()` | Load Bruker raw spectra from a directory |
| `preprocess_maldi_spectra()` | Smooth and baseline-correct spectra |
| `find_peaks()` / `find_peaks_spectra()` | Detect ordinary and shoulder peaks (single / list) |
| `find_peaks_fast()` / `find_peaks_spectra_fast()` | Fast local peak detection (single / list) |
| `filter_peaks()` / `filter_peaks_spectra()` | Filter peaks by intensity, prominence, and strength |
| `find_frequent_mz()` | Find frequent m/z values across a cohort |
| `build_matched_peaks_matrix()` | Assemble a cohort peak intensity matrix |
| `plot_spectrum()` / `plot_spectra()` | Visualize spectra (requires `ggplot2`) |
| `heatmap_spectrum_matrix()` | Heatmap of a spectrum/peak matrix (requires `pheatmap`) |

### Suggested packages

Some functions load optional packages only when used:

- `plot_spectrum()`, `plot_spectra()` → `ggplot2`
- `heatmap_spectrum_matrix()` → `pheatmap`

---

## Author

**Wonseok Oh** ([ORCID: 0009-0002-0687-8466](https://orcid.org/0009-0002-0687-8466))
