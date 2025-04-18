
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scPaX

**scPaX** is an R package designed to compute and visualize the
probability of accessibility (Paₓ) of regulatory regions from
single‑cell ATAC‑seq (scATAC‑seq) datasets. It uses a balls‑into‑bins
model to fit a nonlinear curve and estimate the key parameter *a*, then
produces a summary table and two publication‑ready plots: \* A **scatter
plot** of Paₓ vs. number of fragments per region \* A **histogram of**
log₁₀(Paₓ) distribution

## Installation

Install development tools (if not already installed):

``` r
install.packages(c("devtools", "usethis"))
```

Then install scPaX directly from GitHub:

``` r
devtools::install_github("elesanesc/scPaX")
```

Or, if you are working inside the package directory:

``` r
# From within the scPaX project root
devtools::install()
```

## Usage

``` r
library(scPaX)
library(Seurat)
library(Signac)

# Load your processed Seurat object containing the scATAC assay
seurat_obj <- readRDS("path/to/your_seurat_object.rds")

# Run PaX calculation
results <- calculate_PaX(
  seurat_obj = seurat_obj,
  assay      = "peaks"   # name of the assay slot
)

# Estimated parameter 'a'
results$a
#> [1] 0.002345

# View the top rows of the results table
head(results$data)

# View the scatter plot and histogram
pax_results$scatter_plot
pax_results$histogram
```

## Saving Outputs

After running the function, you can save the output table as a .csv,
.tsv, etc, and the plots as .png, .jpg, .pdf, etc.

``` r
# Save the table to CSV
write.csv(pax_results$data, "PaX_results.csv", row.names = FALSE)

# Save the plots
ggsave("PaX_scatter_plot.png", pax_results$scatter_plot, width = 8, height = 6)
ggsave("PaX_histogram.png", pax_results$histogram, width = 8, height = 6)
```

## Funtion Specifications

``` r
calculate_PaX(seurat_obj = seurat_obj, assay = "peaks")
```

#### Arguments

  - *seurat\_obj*: A Seurat object with a single‑cell ATAC assay.
  - *assay*: Name of the assay slot (default “peaks”).

#### Value

Returns a list with components: \* *model*: the fitted nonlinear model
(nls object) \* *a*: numeric estimate of the parameter *a* \* *data*:
data frame with fragment counts, raw & fitted Paₓ per region \*
*scatter\_plot*: ggplot2 object of Paₓ vs. fragments \* *histogram*:
ggplot2 object of log₁₀(Paₓ) distribution

## Contact

If you have any comments or suggestions about scPaX please raise an
issue or contact us: Elena Sánchez-Escabias:
<elena.sanchez.escabias@cabimer.es> JC Reyes: <jose.reyes@cabimer.es>

## License

This project is licensed under the MIT License. See the LICENSE file for
details.
