# This script plots all treestats returned by get_all_treestats.R

library(gridExtra)
library(grid)
library(ggplot2)

fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
out_fp          = paste0(fp,"plot/sim_treestats/")

# load tree stats data from Fig 4 (large dd interval)
mesh_wd = read.table(paste0(out_fp,'mesh_wd.csv'),sep = ",",header = T)
mesh_ed = read.table(paste0(out_fp,'mesh_ed.csv'),sep = ",",header = T)
mesh_dd = read.table(paste0(out_fp,'mesh_dd.csv'),sep = ",",header = T)

# plot
# match the y-limit for each stat across dd-within, dd-ext, dd-disp
y_limits = lapply(common_stats, function(stat) {
  values = c(mesh_wd[[stat]],
             mesh_dd[[stat]],
             mesh_ed[[stat]])
  range_bound = range(values, na.rm = TRUE)
  # if any of the value is non-negative, set lower bound = 0
  if(range_bound[1] >= 0){
    range_bound[1] = 0
  }
  range_bound
})

names(y_limits) = common_stats

# plot function to create all tree stat plots for a given submodel
make_plots = function(dat, x_col, plot_label, y_limits) {
  #
  label_map <- c(
    w_d      = "w^D",
    e_d      = "e^D",
    d_d_away = 'd^{D*","*dest}'
  )
  lab = label_map[x_col]
  # Select numeric statistic columns, excluding the x variable and other unused variables
  stat_cols = names(dat)[vapply(dat, is.numeric, logical(1))]
  stat_cols = setdiff(stat_cols, c("idx","w_d","e_d","d_d_away","prop_invalid_trees"))
  
  plots <- lapply(stat_cols, function(stat) {
    p = ggplot(dat,aes(x = .data[[x_col]], y = .data[[stat]])) 
    p = p + geom_line(
      color = "skyblue",
      linewidth = 1,
      na.rm = TRUE
    ) +
      labs(
        title = paste(plot_label, stat),
        x = parse(text = lab),
        y = stat
      ) 
    if (x_col == "w_d" || x_col == "d_d_away"){
      p = p + coord_cartesian(ylim = y_limits[[stat]])+ theme_minimal()
    } else {
      p = p + coord_cartesian(ylim = y_limits[[stat]])+ theme_minimal() + scale_x_reverse()
    }
    p
  })
  
  names(plots) = stat_cols
  
  plots
}

# dd-ext
plots_ext = make_plots(dat = mesh_ed,x_col = "e_d",plot_label = "Extinction:",y_limits = y_limits)
# dd-within
plots_within = make_plots(dat = mesh_wd,x_col = "w_d",plot_label = "Speciation:",y_limits = y_limits)
# dd-dispersal
plots_dispersal = make_plots(dat = mesh_dd,x_col = "d_d_away",plot_label = "Dispersal:",y_limits = y_limits)

# find all common stats across all submodels
common_stats = Reduce(
  intersect,
  list(
    names(plots_ext),
    names(plots_within),
    names(plots_dispersal)
  )
)

# combine plots for a given treestat across all submodels, display them in the same window
combined_plots = lapply(common_stats, function(stat) {
  
  gridExtra::arrangeGrob(
    plots_ext[[stat]],
    plots_within[[stat]],
    plots_dispersal[[stat]],
    nrow = 2,
    ncol = 2,
    top = stat
  )
})

names(combined_plots) = common_stats

#save
pdf(
  file = paste0(out_fp,"all_treestats_dd_large.pdf"),
  width = 8,
  height = 6
)

for (stat in common_stats) {
  gridExtra::grid.arrange(combined_plots[[stat]])
}

dev.off()
