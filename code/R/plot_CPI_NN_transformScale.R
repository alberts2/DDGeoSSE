# This script visualizes the CPI plot from NN estimation on each DDGeoSSE parameter
# in transformed scale

### LIBRARY ###
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

### FILE PATH ###
path = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Anolis/final_submodel_sum_no_weight/"
path_in = paste0(path,"submodel_7_mixture_anolis/estimate/")


# ============================================================
# Read and combine the true and estimated parameter values
# ============================================================

true_df = read.csv(paste0(path_in,"sim.test_true.labels_num.csv"))
est_df  = read.csv(paste0(path_in,"sim.test_est.labels_num.csv"))

# Match true and estimated values by simulation index
dat <- inner_join(true_df, est_df, by = "idx")


# ============================================================
# Function for making a CPI plot
# ============================================================

make_accuracy_plot =  function(dat,
                               true_col,
                               est_col,
                               lower_col,
                               upper_col,
                               param_label,
                               n_show = 50) {
  
  # Extract the relevant columns for this parameter
  d = data.frame(true     = dat[[true_col]],
                 estimate = dat[[est_col]],
                 lower    = dat[[lower_col]],
                 upper    = dat[[upper_col]]
                )
  
  # Remove incomplete observations, if any
  d = d[complete.cases(d), ]
  
  # ----------------------------------------------------------
  # Statistics calculated using ALL test examples
  # ----------------------------------------------------------
  
  mae = mean(abs(d$estimate - d$true))
  
  mse = mean((d$estimate - d$true)^2)
  
  rmse = sqrt(mse)
  
  # Whether the 80% CPI contains the true value
  d$covered = (d$lower < d$true) &
    (d$upper > d$true)
  
  coverage = mean(d$covered) * 100
  
  
  # Linear regression using ALL test examples
  fit = lm(estimate ~ true, data = d)
  
  intercept = coef(fit)[1]
  slope     = coef(fit)[2]
  
  # ----------------------------------------------------------
  # First 50 examples for visualization
  # This matches phyddle: these are NOT randomly sampled
  # ----------------------------------------------------------
  
  n_plot = min(n_show, nrow(d))
  
  d_show = d[seq_len(n_plot), ]
  
  
  # ----------------------------------------------------------
  # Axis limits based on the displayed examples
  # ----------------------------------------------------------
  
  min_lim = min(
    c(d_show$true, d_show$lower),
    na.rm = TRUE
  )
  
  max_lim = max(
    c(d_show$true, d_show$upper),
    na.rm = TRUE
  )
  
  range_lim = max_lim - min_lim
  
  min_lim = min_lim - 0.05 * range_lim
  max_lim = max_lim + 0.05 * range_lim
  
  
  # ----------------------------------------------------------
  # Statistics shown in upper-left corner
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
    aes(x = true, y = estimate)
  ) +
    
    # 80% CPI
    geom_errorbar(
      aes(
        ymin = lower,
        ymax = upper,
        color = covered
      ),
      width = 0.08,
      linewidth = 0.4,
      alpha = 0.5
    ) +
    
    # Point estimates
    geom_point(
      aes(color = covered),
      size = 1.2,
      alpha = 0.7
    ) +
    
    # 1:1 line
    geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 0.7
    ) +
    
    # Regression line calculated using ALL test examples
    geom_abline(
      intercept = intercept,
      slope = slope,
      linetype = "dotted",
      linewidth = 0.7
    ) +
    
    # Statistics
    annotate(
      "text",
      x = min_lim,
      y = max_lim,
      label = stat_text,
      hjust = 0,
      vjust = 1,
      size = 3.9,
      lineheight = 0.9
    ) +
    
    # Large coverage number in bottom-right
    annotate(
      "text",
      x = max_lim,
      y = min_lim,
      label = sprintf("%.1f%%", coverage),
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
      xlim = c(min_lim, max_lim),
      ylim = c(min_lim, max_lim)
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
        size = 20
      ),
      
      axis.text = element_text(
        size = 14
      ),
      
      panel.grid = element_blank(),
      # remove separate axis lines
      axis.line = element_blank(),
      
      # full rectangular box
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      )
    )
  
  return(p)
}


# ============================================================
# Make the seven panels
# ============================================================

p_w = make_accuracy_plot(
  dat,
  true_col  = "log_base_w",
  est_col   = "log_base_w_value",
  lower_col = "log_base_w_lower",
  upper_col = "log_base_w_upper",
  param_label = expression(log[e](rho[w]))
)


p_e = make_accuracy_plot(
  dat,
  true_col  = "log_base_e",
  est_col   = "log_base_e_value",
  lower_col = "log_base_e_lower",
  upper_col = "log_base_e_upper",
  param_label = expression(log[e](rho[e]))
)


p_d = make_accuracy_plot(
  dat,
  true_col  = "log_base_d",
  est_col   = "log_base_d_value",
  lower_col = "log_base_d_lower",
  upper_col = "log_base_d_upper",
  param_label = expression(log[e](rho[d]))
)


p_b = make_accuracy_plot(
  dat,
  true_col  = "log_base_b",
  est_col   = "log_base_b_value",
  lower_col = "log_base_b_lower",
  upper_col = "log_base_b_upper",
  param_label = expression(log[e](rho[b]))
)


p_eD = make_accuracy_plot(
  dat,
  true_col  = "exp_scalar_e",
  est_col   = "exp_scalar_e_value",
  lower_col = "exp_scalar_e_lower",
  upper_col = "exp_scalar_e_upper",
  param_label = expression(exp(e^D))
)


p_wD = make_accuracy_plot(
  dat,
  true_col  = "exp_scalar_w",
  est_col   = "exp_scalar_w_value",
  lower_col = "exp_scalar_w_lower",
  upper_col = "exp_scalar_w_upper",
  param_label = expression(exp(w^D))
)


p_dD = make_accuracy_plot(
  dat,
  true_col  = "exp_scalar_d_j",
  est_col   = "exp_scalar_d_j_value",
  lower_col = "exp_scalar_d_j_lower",
  upper_col = "exp_scalar_d_j_upper",
  param_label = expression(exp(d^{D*","*dest}))
)


# ============================================================
# Arrange panels
# ============================================================

top_row = (
  p_w + p_e + p_d + p_b
) +
  plot_layout(ncol = 4)


bottom_row = (
  plot_spacer() +
    p_eD +
    p_wD +
    p_dD +
    plot_spacer()
) +
  plot_layout(
    widths = c(0.5, 1, 1, 1, 0.5)
  )


fig5 = top_row / bottom_row

final_fig = ggdraw(fig5) +
  draw_label(
    "CNN prediction",
    x = 0.015,
    y = 0.5,
    angle = 90,
    size = 18
  ) +
  draw_label(
    "True value",
    x = 0.5,
    y = 0.015,
    size = 18
  )

final_fig