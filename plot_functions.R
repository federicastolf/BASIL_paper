
# Function to create correlation heatmaps 
plot_correlation_heatmap = function(data, title, cluster = FALSE, row_order = NULL, 
                                    col_order = NULL) {
  pheatmap(mat = if (!is.null(row_order) && !is.null(col_order)) 
    data[row_order, col_order] else data,
    treeheight_row = 0, treeheight_col = 0, cluster_rows = cluster, 
    cluster_cols = cluster, border_color = NA, legend = FALSE, show_colnames = FALSE,
    show_rownames = FALSE, breaks = breaks, color = cols, fontsize = 12, main = title)
}

# Function to create scatter plot with identity line
plot_correlation_scatter = function(data, title, lim_ax, point_color = "#1170aa") {
  ggplot(data, aes(x = Observed, y = Predicted)) +
    geom_point(color = point_color, shape = 1) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1) +
    labs(title = title, x = "Observed", y = "Predicted") +
    xlim(lim_ax) + ylim(lim_ax) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_text(size = 16))
}
