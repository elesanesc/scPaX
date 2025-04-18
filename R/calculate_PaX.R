#' scPaX: Calculate Paₓ (Probability of Accessibility) for scATAC-seq data
#'
#' This function calculates Paₓ values from a scATAC-seq Seurat object, fits a nonlinear model
#' to estimate parameter `a`, and generates summary plots and a result table.
#'
#' @param seurat_obj A Seurat object with a scATAC assay (default assay="peaks").
#' @param assay The name of the assay containing the peak count matrix. Default: "peaks".
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted nonlinear model object.}
#'   \item{a}{Estimated value of the parameter `a`.}
#'   \item{data}{A data frame with peak-level statistics and Paₓ estimates.}
#'   \item{scatter_plot}{A ggplot object for the Paₓ vs fragments scatter plot.}
#'   \item{histogram}{A ggplot object for the log10(Paₓ) histogram.}
#' }
#'
#' @importFrom stats nls coef predict
#' @importFrom ggplot2 ggplot aes geom_point labs ylim annotate theme_classic geom_histogram
#' @importFrom utils write.csv
#'
#' @export
calculate_PaX <- function(seurat_obj, assay = "peaks") {

  message("Extracting count matrix from assay: ", assay)

  # Extract the count matrix from the specified assay
  atac_counts <- as.matrix(seurat_obj@assays[[assay]]@counts)

  # Filter peaks with at least one fragment in at least one cell
  binary_matrix <- atac_counts
  binary_matrix[binary_matrix > 1] <- 1
  atac_filtered <- atac_counts[rowSums(binary_matrix) > 0, ]

  # Compute peak-level statistics
  message("Computing fragment and cell counts per peak...")
  cells_with_at_least_1_frag <- rowSums(atac_filtered > 0)
  cells_with_at_least_2_frags <- rowSums(atac_filtered >= 2)
  total_fragments <- rowSums(atac_filtered)

  peak_summary <- data.frame(
    peak_id = rownames(atac_filtered),
    cells_with_at_least_1_frag = cells_with_at_least_1_frag,
    cells_with_at_least_2_frags = cells_with_at_least_2_frags,
    total_fragments = total_fragments,
    fragments_per_cell = total_fragments / cells_with_at_least_1_frag
  )

  # Model setup
  message("Fitting non-linear model to estimate 'a' in Pa\u2093 = 1 - exp(-m\u2093/n) - (a * m\u2093)/n ...")
  num_cells <- ncol(atac_filtered)

  model_data <- data.frame(
    Pa_x = peak_summary$cells_with_at_least_1_frag / num_cells,
    m_x = peak_summary$total_fragments
  )

  # Fit non-linear least squares model
  model_fit <- nls(Pa_x ~ 1 - exp(-m_x / num_cells) - (a * m_x) / num_cells,
                   data = model_data,
                   start = list(a = 1))

  a_estimate <- coef(model_fit)["a"]
  Pa_x_fitted <- predict(model_fit)

  message("Model fit complete. Estimated a = ", round(a_estimate, 6))

  # Add predictions to summary table
  peak_summary$Pa_x_fitted <- Pa_x_fitted

  message("Generating plots...")
  summary_ordered <- peak_summary[order(peak_summary$total_fragments), ]
  pa_min <- round(min(summary_ordered$Pa_x_fitted, na.rm = TRUE), 6)
  pa_max <- round(max(summary_ordered$Pa_x_fitted, na.rm = TRUE), 6)

  # Scatter plot
  scatter_plot <- ggplot(summary_ordered, aes(x = total_fragments, y = Pa_x_fitted)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    labs(
      x = "Number of fragments of AR\u2093",
      y = expression(italic(Pa)[x]),
      title = "Probability of accessibility for AR\u2093"
    ) +
    ylim(0, 1) +
    annotate("text", x = 1, y = 1,
             label = paste0(pa_min, " <= Pa\u2093 <= ", pa_max),
             hjust = 0, vjust = 1, size = 5) +
    theme_classic(base_size = 14)

  # Histogram
  filtered_summary <- subset(peak_summary, Pa_x_fitted > 0)
  histogram <- ggplot(filtered_summary, aes(x = log10(Pa_x_fitted))) +
    geom_histogram(binwidth = binwidth, fill = "steelblue", color = "black") +
    labs(
      x = expression(Log[10](Pa[x])),
      y = "Number of AR\u2093s",
      title = "Distribution of Pa\u2093 frequencies"
    ) +
    theme_classic(base_size = 14)

  message("Done!")

  # Return results
  return(list(
    model = model_fit,
    a = a_estimate,
    data = peak_summary,
    scatter_plot = scatter_plot,
    histogram = histogram
  ))
}
#' @examples
#' pax_results <- calculate_PaX(seurat_obj = seurat_obj, output_prefix = "/my_path/")
#'
#' # Save the table to CSV
#' write.csv(pax_results$data, "PaX_results.csv", row.names = FALSE)
#'
#' # Save the plots
#' ggsave("PaX_scatter_plot.png", pax_results$scatter_plot, width = 8, height = 6)
#' ggsave("PaX_histogram.png", pax_results$histogram, width = 8, height = 6)
