library(ggplot2)
library(gridExtra)
library(grid)

fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
out_fp          = paste0(fp,"plot/sim_treestats/")

which_dat  = "small" # large (fig.4) or small (supplementary figure)
if (which_dat == "large"){
  # load tree stats data from Fig 4 (large dd interval)
  mesh_wd = read.table(paste0(out_fp,'mesh_wd_big.csv'),sep = ",",header = T)
  mesh_ed = read.table(paste0(out_fp,'mesh_ed_big.csv'),sep = ",",header = T)
  mesh_dd = read.table(paste0(out_fp,'mesh_dd_big.csv'),sep = ",",header = T)
} else if (which_dat == "small"){
  # load tree stats data from supplement (small dd interval)
  mesh_wd = read.table(paste0(out_fp,'ladder_wd_small.csv'),sep = ",",header = T)
  mesh_ed = read.table(paste0(out_fp,'ladder_ed_small.csv'),sep = ",",header = T)
  mesh_dd = read.table(paste0(out_fp,'ladder_dd_small.csv'),sep = ",",header = T)
}


# stats of interest 
stats = c(
            "colless_corr",
            "gamma",
            "mean_extant_taxa",
            "mean_range_size",
            "avg_ladder"
          )

# match the ylimits in the same column
y_limits = list()

for (stat in stats) {
  
  values = c(
    mesh_wd[[stat]],
    mesh_ed[[stat]],
    mesh_dd[[stat]]
  )
  
  y_limits[[stat]] <- range(values, na.rm = TRUE)
}


#---------------------------------------
# Function for making an individual plot
#---------------------------------------
make_panel = function(
    dat,
    x_col,
    y_col,
    y_lim = NULL,
    line_colour,
    x_label = NULL,
    title = NULL,
    reverse_x = FALSE,
    log_y = FALSE
) {
  
  p = ggplot(
    dat,
    aes(x = .data[[x_col]], y = .data[[y_col]])
  ) +
    geom_line(
      colour = line_colour,
      linewidth = 1
    ) +
    labs(
      title = title,
      x = x_label,
      y = NULL
    ) +
    theme_minimal(base_size = 12) + 
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text = element_text(size = 15),
      plot.margin = margin(5, 8, 5, 8)
    )

  if (which_dat == "small") {
    
    small_x_breaks = seq(
      min(dat[[x_col]], na.rm = TRUE),
      max(dat[[x_col]], na.rm = TRUE),
      length.out = 3
    )
    
    if (reverse_x) {
      p = p + scale_x_reverse(
        breaks = small_x_breaks,
        expand = expansion(mult = c(0.05, 0.10))
      )
    } else {
      p = p + scale_x_continuous(
        breaks = small_x_breaks,
        expand = expansion(mult = c(0.05, 0.10))
      )
    }
    
  } else if (reverse_x) {
    p = p + scale_x_reverse()
  }
  
  if (log_y) {
    p = p +
      scale_y_log10(
        limits = y_lim,
        breaks = c(1, 10, 100, 1000)
      )
  } else if (!is.null(y_lim)) {
    p = p + coord_cartesian(ylim = y_lim)
  }
  
  p
}

# DD-within row

p_wd_ladder = make_panel(
  mesh_wd,
  x_col = "w_d",
  y_col = "avg_ladder",
  y_lim  = y_limits[["avg_ladder"]],
  line_colour = "orange",
  title = "Mean ladder",
  reverse_x = TRUE
)

p_wd_colless = make_panel(
  mesh_wd,
  x_col = "w_d",
  y_col = "colless_corr",
  y_lim  = y_limits[["colless_corr"]],
  line_colour = "skyblue",
  title = "Mean corrected Colless",
  reverse_x = TRUE
)

p_wd_gamma = make_panel(
  mesh_wd,
  x_col = "w_d",
  y_col = "gamma",
  y_lim  = y_limits[["gamma"]],
  line_colour = "lightgreen",
  title = expression(bold("Mean") ~ gamma),
  reverse_x = TRUE
)

p_wd_taxa = make_panel(
  mesh_wd,
  x_col = "w_d",
  y_col = "mean_extant_taxa",
  y_lim  = y_limits[["mean_extant_taxa"]],
  line_colour = "black",
  title = "Mean extant taxa",
  log_y = TRUE,
  reverse_x = TRUE
)

p_wd_range = make_panel(
  mesh_wd,
  x_col = "w_d",
  y_col = "mean_range_size",
  y_lim  = y_limits[["mean_range_size"]],
  line_colour = "red",
  title = "Mean range size",
  reverse_x = TRUE
)

# DD-ext row

p_ed_ladder = make_panel(
  mesh_ed,
  x_col = "e_d",
  y_col = "avg_ladder",
  y_lim  = y_limits[["avg_ladder"]],
  line_colour = "orange",
  reverse_x = FALSE
)

p_ed_colless = make_panel(
  mesh_ed,
  x_col = "e_d",
  y_col = "colless_corr",
  y_lim = y_limits[["colless_corr"]],
  line_colour = "skyblue",
  x_label = expression(e^D),
  reverse_x = FALSE
)

p_ed_gamma = make_panel(
  mesh_ed,
  x_col = "e_d",
  y_col = "gamma",
  y_lim = y_limits[["gamma"]],
  line_colour = "lightgreen",
  x_label = expression(e^D),
  reverse_x = FALSE
)

p_ed_taxa = make_panel(
  mesh_ed,
  x_col = "e_d",
  y_col = "mean_extant_taxa",
  y_lim = y_limits[["mean_extant_taxa"]],
  line_colour = "black",
  x_label = expression(e^D),
  reverse_x = FALSE,
  log_y = TRUE
)

p_ed_range = make_panel(
  mesh_ed,
  x_col = "e_d",
  y_col = "mean_range_size",
  y_lim = y_limits[["mean_range_size"]],
  line_colour = "red",
  x_label = expression(e^D),
  reverse_x = FALSE
)

# DD-dispersal row

p_dd_ladder = make_panel(
  mesh_dd,
  x_col = "d_d_away",
  y_col = "avg_ladder",
  y_lim  = y_limits[["avg_ladder"]],
  line_colour = "orange",
  reverse_x = TRUE
)

p_dd_colless = make_panel(
  mesh_dd,
  x_col = "d_d_away",
  y_col = "colless_corr",
  y_lim = y_limits[["colless_corr"]],
  line_colour = "skyblue",
  x_label = expression(d^{D*","*dest}),
  reverse_x = TRUE
)

p_dd_gamma = make_panel(
  mesh_dd,
  x_col = "d_d_away",
  y_col = "gamma",
  y_lim = y_limits[["gamma"]],
  line_colour = "lightgreen",
  x_label = expression(d^{D*","*dest}),
  reverse_x = TRUE
)

p_dd_taxa = make_panel(
  mesh_dd,
  x_col = "d_d_away",
  y_col = "mean_extant_taxa",
  y_lim = y_limits[["mean_extant_taxa"]],
  line_colour = "black",
  x_label = expression(d^{D*","*dest}),
  log_y = TRUE,
  reverse_x = TRUE
)

p_dd_range = make_panel(
  mesh_dd,
  x_col = "d_d_away",
  y_col = "mean_range_size",
  y_lim = y_limits[["mean_range_size"]],
  line_colour = "red",
  x_label = expression(d^{D*","*dest}),
  reverse_x = TRUE
)

# Remove individual x-axis titles from every panel
p_wd_ladder  = p_wd_ladder  + labs(x = NULL)
p_wd_colless = p_wd_colless + labs(x = NULL)
p_wd_gamma   = p_wd_gamma   + labs(x = NULL)
p_wd_taxa    = p_wd_taxa    + labs(x = NULL)
p_wd_range   = p_wd_range   + labs(x = NULL)

p_ed_ladder  = p_ed_ladder  + labs(x = NULL)
p_ed_colless = p_ed_colless + labs(x = NULL)
p_ed_gamma   = p_ed_gamma   + labs(x = NULL)
p_ed_taxa    = p_ed_taxa    + labs(x = NULL)
p_ed_range   = p_ed_range   + labs(x = NULL)

p_dd_ladder  = p_dd_ladder  + labs(x = NULL)
p_dd_colless = p_dd_colless + labs(x = NULL)
p_dd_gamma   = p_dd_gamma   + labs(x = NULL)
p_dd_taxa    = p_dd_taxa    + labs(x = NULL)
p_dd_range   = p_dd_range   + labs(x = NULL)

# DD-within row
row_wd_plots <- arrangeGrob(
  p_wd_taxa,
  p_wd_range,
  p_wd_colless,
  p_wd_ladder,
  p_wd_gamma,
  nrow = 1
)
# plot for row wd
row_wd = arrangeGrob(
  row_wd_plots,
  textGrob(
    expression(w^D),
    gp = gpar(fontsize = 16)
  ),
  ncol = 1,
  heights = c(1, 0.15)
)

# DD-extinction row
row_ed_plots = arrangeGrob(
  p_ed_taxa,
  p_ed_range,
  p_ed_colless,
  p_ed_ladder,
  p_ed_gamma,
  nrow = 1
)
# plot for row ed
row_ed <- arrangeGrob(
  row_ed_plots,
  textGrob(
    expression(e^D),
    gp = gpar(fontsize = 16)
  ),
  ncol = 1,
  heights = c(1, 0.15)
)

# DD-dispersal row
row_dd_plots = arrangeGrob(
  p_dd_taxa,
  p_dd_range,
  p_dd_colless,
  p_dd_ladder,
  p_dd_gamma,
  nrow = 1
)
# plot for row dd
row_dd = arrangeGrob(
  row_dd_plots,
  textGrob(
    expression(d^{D*","*dest}),
    gp = gpar(fontsize = 16)
  ),
  ncol = 1,
  heights = c(1, 0.15)
)

# label y-axis
label_wd = textGrob(
  "DD-within",
  rot = 90,
  gp = gpar(
    fontsize = 16,
    fontface = "bold"
  )
)

label_ed = textGrob(
  "DD-extinction",
  rot = 90,
  gp = gpar(
    fontsize = 16,
    fontface = "bold"
  )
)

label_dd = textGrob(
  "DD-dispersal",
  rot = 90,
  gp = gpar(
    fontsize = 16,
    fontface = "bold"
  )
)

#
complete_row_wd = arrangeGrob(
  label_wd,
  row_wd,
  ncol = 2,
  widths = c(0.06, 1)
)

complete_row_ed <- arrangeGrob(
  label_ed,
  row_ed,
  ncol = 2,
  widths = c(0.06, 1)
)

complete_row_dd <- arrangeGrob(
  label_dd,
  row_dd,
  ncol = 2,
  widths = c(0.06, 1)
)

combined_plot <- arrangeGrob(
  complete_row_wd,
  complete_row_ed,
  complete_row_dd,
  ncol = 1,
  heights = c(1, 1, 1)
)

# Save
grid.newpage()
grid.draw(combined_plot)

if (which_dat == "large"){
  ggsave(
    filename = paste0(
      out_fp,
      "treestats_dd_all_large.pdf"
    ),
    plot = combined_plot,
    width = 15,
    height = 9
  )
} else if (which_dat == "small"){
  ggsave(
    filename = paste0(
      out_fp,
      "treestats_dd_all_small.pdf"
    ),
    plot = combined_plot,
    width = 15,
    height = 9
  )
}

