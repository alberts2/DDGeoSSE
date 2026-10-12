# This script visualizes the CPI plot from NN estimation on each DDGeoSSE parameter
# in original scale (not transformed scale)

### LIBRARY ###
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(grid)

### FILE PATH ###
path     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Anolis/final_submodel_sum_no_weight/"

out_base = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/plot/deep_estimate/"


# ============================================================
# Function for making a CPI plot
# ============================================================

make_accuracy_plot = function(dat,
                               true_col,
                               est_col,
                               lower_col,
                               upper_col,
                               param_label,
                               transform = "none",
                               n_show = 50,
                               fixed_limits = NULL,
                               submodel = NA,
                               param_name = "") {
  
  ### Extract corresponding columns from dataset ###
  
  d = data.frame(original_row = seq_len(nrow(dat)),
                 true         = dat[[true_col]],
                 estimate     = dat[[est_col]],
                 lower        = dat[[lower_col]],
                 upper        = dat[[upper_col]]
                 )
  
  # Remove pre-existing missing values
  d = d[complete.cases(d), ]
  
  
  # ----------------------------------------------------------
  # Back-transform to ORIGINAL DDGeoSSE parameter scale
  # ----------------------------------------------------------
  
  if (transform == "exp") {
    
    # Network was trained on log(parameter) e.g., base rate parameters
    # log(rho) = rho
    
    d$true     = exp(d$true)
    d$estimate = exp(d$estimate)
    d$lower    = exp(d$lower)
    d$upper    = exp(d$upper)
    
    
  } else if (transform == "log") {
    
    # Network was trained on exp(D) e.g., DD effect parameters
    # exp(D) = D
    #
    # log() requires all quantities to be strictly positive.
    
    invalid = (d$true <= 0 |
                d$estimate <= 0 |
                d$lower <= 0 |
                d$upper <= 0
                )
    
    n_invalid = sum(invalid, na.rm = TRUE)
  
    
    # Report any observations that cannot be back-transformed
    
    if (n_invalid > 0) {
      
      cat(
        "\n============================================\n",
        "Non-positive value before log()\n",
        "Submodel: ", submodel, "\n",
        "Parameter: ", param_name, "\n",
        "Excluded: ", n_invalid, " of ", nrow(d), "\n",
        "============================================\n",
        sep = ""
      )
      
      print(
        data.frame(
          d[invalid, , drop = FALSE]
        )
      )
    }
    
    
    # Remove observations that cannot be mapped back to D scale
    
    d = d[!invalid, , drop = FALSE]
    
    
    # Back-transform to original DDGeoSSE scale
    
    d$true     = log(d$true)
    d$estimate = log(d$estimate)
    d$lower    = log(d$lower)
    d$upper    = log(d$upper)
  }
  
  
  # ----------------------------------------------------------
  # Statistics on ORIGINAL scale using ALL valid test examples
  # ----------------------------------------------------------
  
  mae = mean(
    abs(d$estimate - d$true)
  )
  
  mse = mean(
    (d$estimate - d$true)^2
  )
  
  rmse = sqrt(mse)
  
  
  # Coverage
  
  d$covered = (d$lower < d$true & d$upper > d$true)
  
  coverage = mean(d$covered) * 100
  
  
  # Regression
  
  fit = lm(
    estimate ~ true,
    data = d
  )
  
  intercept = coef(fit)[1]
  slope     = coef(fit)[2]
  
  
  # ----------------------------------------------------------
  # Same test examples as phyddle:
  # original first 50 rows only
  # ----------------------------------------------------------
  
  d_show = d[
    d$original_row <= n_show,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Axis limits
  # ----------------------------------------------------------
  
  if (is.null(fixed_limits)) {
    
    min_lim = min(
      c(d_show$true,
        d_show$lower
      ),
      na.rm = TRUE
    )
    
    max_lim = max(
      c(d_show$true,
        d_show$upper
      ),
      na.rm = TRUE
    )
    
    range_lim = max_lim - min_lim
    
    min_lim = min_lim - 0.05 * range_lim
    max_lim = max_lim + 0.05 * range_lim
    
  } else {
    
    min_lim = fixed_limits[1]
    max_lim = fixed_limits[2]
    
    range_lim = max_lim - min_lim
  }
  
  
  # Error-bar cap width proportional to x-axis range
  
  cap_width = 0.015 * range_lim
  
  
  # ----------------------------------------------------------
  # Statistics label
  # ----------------------------------------------------------
  
  stat_text = paste0(
    "MAE: ", sprintf("%.2E", mae), "\n",
    "MSE: ", sprintf("%.2E", mse), "\n",
    "RMSE: ", sprintf("%.2E", rmse), "\n",
    "Intercept: ", sprintf("%.2E", intercept), "\n",
    "Slope: ", sprintf("%.2E", slope), "\n",
    "Coverage: ", sprintf("%.1f", coverage), "%\n",
    "Coverage target: 80.0%"
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p = ggplot(
    d_show,
    aes(
      x = true,
      y = estimate
    )
  ) +
    
    # CPI
    
    geom_errorbar(
      aes(
        ymin = lower,
        ymax = upper,
        color = covered
      ),
      width = cap_width,
      linewidth = 0.4,
      alpha = 0.5
    ) +
    
    # Point estimate
    
    geom_point(
      aes(
        color = covered
      ),
      size = 1.2,
      alpha = 0.7
    ) +
    
    # 1:1 line
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 0.7
    ) +
    
    # Regression line
    
    geom_abline(
      intercept = intercept,
      slope = slope,
      linetype = "dotted",
      linewidth = 0.7
    ) +
    
    # Statistics
    
    annotate(
      "label",
      x = min_lim + 0.03 * range_lim,
      y = max_lim - 0.03 * range_lim,
      label = stat_text,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      lineheight = 0.9,
      fill = "white",
      label.size = NA,
      label.padding = unit(
        0.05,
        "lines"
      )
    ) +
    
    # Large coverage value
    
    annotate(
      "text",
      x = max_lim,
      y = min_lim,
      label = sprintf(
        "%.1f%%",
        coverage
      ),
      hjust = 1,
      vjust = 0,
      size = 7
    ) +
    
    scale_color_manual(
      values = c(
        "TRUE"  = "purple",
        "FALSE" = "grey80"
      ),
      guide = "none"
    ) +
    
    coord_fixed(
      xlim = c(
        min_lim,
        max_lim
      ),
      ylim = c(
        min_lim,
        max_lim
      )
    ) +
    
    labs(
      title = param_label,
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal() +
    
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        size = 20,
        margin = margin(
          b = 5
        )
      ),
      
      axis.text = element_text(
        size = 14
      ),
      
      panel.grid = element_blank(),
      
      axis.line = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      )
    )
  
  return(p)
}


# ============================================================
# Loop over submodels
# ============================================================

for (submodel in 0:7) {
  
  cat(
    "\n\n############################################\n",
    "Making figure for submodel ",
    submodel,
    "\n############################################\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Input/output directories
  # ----------------------------------------------------------
  
  path_in <- paste0(
    path,
    "submodel_",
    submodel,
    "_mixture_anolis/estimate/"
  )
  
  out_path <- paste0(
    out_base,
    "submodel_",
    submodel,
    "/"
  )
  
  dir.create(
    out_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Read data
  # ----------------------------------------------------------
  
  true_df = read.csv(
    paste0(
      path_in,
      "sim.test_true.labels_num.csv"
    )
  )
  
  est_df = read.csv(
    paste0(
      path_in,
      "sim.test_est.labels_num.csv"
    )
  )
  
  
  # Match true and estimated values by simulation index
  
  dat = inner_join(
    true_df,
    est_df,
    by = "idx"
  )
  
  
  # ==========================================================
  # Base-rate parameters
  # ==========================================================
  
  base_limits = c(0,1.5)
  
  # rho_w
  
  p_w = make_accuracy_plot(
    dat,
    true_col  = "log_base_w",
    est_col   = "log_base_w_value",
    lower_col = "log_base_w_lower",
    upper_col = "log_base_w_upper",
    param_label = expression(rho[w]),
    transform = "exp",
    fixed_limits = base_limits,
    submodel = submodel,
    param_name = "rho_w"
  )
  
  # rho_e
  
  p_e = make_accuracy_plot(
    dat,
    true_col  = "log_base_e",
    est_col   = "log_base_e_value",
    lower_col = "log_base_e_lower",
    upper_col = "log_base_e_upper",
    param_label = expression(rho[e]),
    transform = "exp",
    fixed_limits = base_limits,
    submodel = submodel,
    param_name = "rho_e"
  )
  
  # rho_d
  
  p_d = make_accuracy_plot(
    dat,
    true_col  = "log_base_d",
    est_col   = "log_base_d_value",
    lower_col = "log_base_d_lower",
    upper_col = "log_base_d_upper",
    param_label = expression(rho[d]),
    transform = "exp",
    fixed_limits = base_limits,
    submodel = submodel,
    param_name = "rho_d"
  )
  
  # rho_b
  
  p_b = make_accuracy_plot(
    dat,
    true_col  = "log_base_b",
    est_col   = "log_base_b_value",
    lower_col = "log_base_b_lower",
    upper_col = "log_base_b_upper",
    param_label = expression(rho[b]),
    transform = "exp",
    fixed_limits = base_limits,
    submodel = submodel,
    param_name = "rho_b"
  )
  
  # ==========================================================
  # DD parameters
  # ==========================================================
  
  dd_plots <- list()
  
  # ----------------------------------------------------------
  # e^D
  # ----------------------------------------------------------
  
  if ("exp_scalar_e" %in% names(dat)) {
    
    dd_plots[["eD"]] = make_accuracy_plot(
      dat,
      true_col  = "exp_scalar_e",
      est_col   = "exp_scalar_e_value",
      lower_col = "exp_scalar_e_lower",
      upper_col = "exp_scalar_e_upper",
      param_label = expression(e^D),
      transform = "log",
      submodel = submodel,
      param_name = "eD"
    )
  }
  
  
  # ----------------------------------------------------------
  # w^D
  # ----------------------------------------------------------
  
  if ("exp_scalar_w" %in% names(dat)) {
    
    dd_plots[["wD"]] = make_accuracy_plot(
      dat,
      true_col  = "exp_scalar_w",
      est_col   = "exp_scalar_w_value",
      lower_col = "exp_scalar_w_lower",
      upper_col = "exp_scalar_w_upper",
      param_label = expression(w^D),
      transform = "log",
      submodel = submodel,
      param_name = "wD"
    )
  }
  
  
  # ----------------------------------------------------------
  # d^{D,dest}
  # ----------------------------------------------------------
  
  if ("exp_scalar_d_j" %in% names(dat)) {
    
    dd_plots[["dD"]] = make_accuracy_plot(
      dat,
      true_col  = "exp_scalar_d_j",
      est_col   = "exp_scalar_d_j_value",
      lower_col = "exp_scalar_d_j_lower",
      upper_col = "exp_scalar_d_j_upper",
      param_label = expression(
        d^{D*","*dest}
      ),
      transform = "log",
      submodel = submodel,
      param_name = "dD,dest"
    )
  }
  
  
  # ==========================================================
  # Arrange panels
  # ==========================================================
  
  top_row = (
    p_w +
      p_e +
      p_d +
      p_b
  ) +
    plot_layout(
      ncol = 4
    )
  
  
  if (length(dd_plots) > 0) {
    
    bottom_row = wrap_plots(
      dd_plots,
      nrow = 1
    )
    
    fig <- top_row / bottom_row
    
  } else {
    
    fig = top_row
  }
  
  
  # ==========================================================
  # Global axis labels
  # ==========================================================
  
  final_fig = ggdraw() +
    
    draw_plot(
      fig,
      x = 0.04,
      y = 0.07,
      width = 0.95,
      height = 0.91
    ) +
    
    draw_label(
      "CNN prediction",
      x = 0.012,
      y = 0.53,
      angle = 90,
      size = 18
    ) +
    
    draw_label(
      "True value",
      x = 0.52,
      y = 0.025,
      size = 18
    )
  
  
  # ==========================================================
  # Save
  # ==========================================================
  
  ggsave(
    paste0(
      out_path,
      "Figure_submodel_",
      submodel,
      ".pdf"
    ),
    final_fig,
    width = 14,
    height = 9,
    units = "in"
  )
}