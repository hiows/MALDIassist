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

## Example

The figure below shows the core workflow on real Bruker MALDI-TOF data: a raw spectrum (gray) is smoothed and baseline-corrected (red), then peaks are detected and filtered (blue).

![MALDIassist workflow example](man/figures/README-example.png)

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hiows/MALDIassist")
```

---

## Quick start

The pipeline follows five steps: **load → preprocess → peak picking → peak filtering → visualization**.

```r
library(MALDIassist)
```

### 1. Load

Point `load_maldi_spectra()` at a directory of Bruker flex data files. It returns a named list of raw spectra.

```r
raw_spectra <- load_maldi_spectra(spectra_dir = "data/")
```

### 2. Preprocess

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

### 3. Peak picking

`find_peaks_spectra()` detects peaks (including shoulder peaks) using a Gaussian KDE approach.

```r
peaks_list <- find_peaks_spectra(
  spectra                    = preprocessed_spectra,
  bw                         = 1,      # KDE bandwidth
  hws_peaks                  = 10,
  weight_type                = "raw",
  cutoff_kappa_peak_strength = 0.3,
  n_cores                    = 4
)
```

### 4. Peak filtering

`filter_peaks_spectra()` removes low-quality peaks by intensity, prominence, and strength cutoffs.

```r
filtered_peaks_list <- filter_peaks_spectra(
  spectra                = preprocessed_spectra,
  peaks_list             = peaks_list,
  cutoff_peak_intensity  = 100,
  cutoff_peak_prominence = 100,
  cutoff_peak_strength   = 0.5,
  normalization_type     = "raw"
)
```

### 5. Visualization

Overlay a raw spectrum (gray), its preprocessed version (red), and the filtered peaks (blue) — this reproduces the figure at the top of this README:

```r
example_range <- c(12000, 15000)
i    <- 1
spec <- raw_spectra[[i]]
pp   <- preprocessed_spectra[[i]]
fp   <- filtered_peaks_list[[i]]

in_range <- function(df) subset(df, df[, 1] >= example_range[1] & df[, 1] <= example_range[2])

plot(
  in_range(spec),
  type = "l", col = "gray70", lwd = 1.5,
  ylim = c(0, max(in_range(spec)[, 2])),
  xlab = "m/z", ylab = "Intensity", frame = FALSE
)
lines(in_range(pp), col = "red", lwd = 2)

fp_r <- fp[fp$mz >= example_range[1] & fp$mz <= example_range[2], ]
segments(fp_r$mz, 0, fp_r$mz, fp_r$intensity, col = "blue", lwd = 1.5)
```

You can also overlay all spectra with `ggplot2`:

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

## How to cite

If you use **MALDIassist** in your research, please cite it. From R you can run:

```r
citation("MALDIassist")
```

A BibTeX entry:

```bibtex
@Manual{maldiassist,
  title  = {MALDIassist: Mathematical Utilities for MALDI-TOF Mass Spectrometry},
  author = {Wonseok Oh},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/hiows/MALDIassist}
}
```

> A DOI will be added here once the release is archived on Zenodo.

## License

MIT © 2026 Wonseok Oh. See [LICENSE.md](LICENSE.md) for details.
