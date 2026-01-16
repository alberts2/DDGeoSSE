library(RevGadgets)
library(deeptime)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)

# filesystem
dir = "/Users/mlandis/projects/DDGeoSSE"

clade_idx = 1
dat_dir = c(
    paste0(dir, "/data/empiric_viburnum_rescaled/oreino5.1"),
    paste0(dir, "/data/empiric_anolis/use_this/anolis.5")
)[clade_idx]
plot_fn = c(
    "/Users/mlandis/projects/DDGeoSSE/plot/tree_trait_figure/oreino5.1.tree_range.pdf",
    "/Users/mlandis/projects/DDGeoSSE/plot/tree_trait_figure/anolis.tree_range.pdf"
)        

# read and format tree
tree_fn = paste0(dat_dir, ".tre")
phy = read.tree(tree_fn)
phy$tip.label = gsub("^'|'$", "", phy$tip.label)
phy$edge.length = phy$edge.length / max(branching.times(phy)) # rescale tree to height 1
num_taxa = length(phy$tip.label)

# read and format tip data
dat_fn = paste0(dat_dir, ".dat.csv")
tip_dat = read.csv(dat_fn)
colnames(tip_dat)[2:6] = LETTERS[1:5] 
reg_lbl = LETTERS[1:5]
tip_dat = read.csv(dat_fn, sep=",", header=T, row.names=1, colClasses=c("character"))
for (i in 1:nrow(tip_dat)) {
    for (j in 1:ncol(tip_dat)) {
        if (tip_dat[i,j] == 0) {
            tip_dat[i,j] = NA
        }
    }
}
for (i in 1:ncol(tip_dat)) {
    tip_dat[ !is.na(tip_dat[,i]), i ] = reg_lbl[i]
    tip_dat[,i] = factor(tip_dat[,i])
} 
colnames(tip_dat)=reg_lbl
tip_dat = data.frame(tip_dat)


# plot settings
tip_fsize = c(4,2)
my_colors = c("firebrick3","gold","deepskyblue","darkorchid","forestgreen") #,"blue","cyan","magenta","orange")
my_labels = list(
    c("A: E/C/W Mexico","B: Oaxaca + Chiapas + Upper CA","Lower CA","D: Northern Andes","E: Central Andes"),
    c("A: Lesser Antilles","B: Puerto Rico","C: Cuba + Caymans","D: Hispaniola","E: Jamaica")
)


# make plot
p2 <- ggtree(phy)

p2 = p2 + geom_tiplab( align=TRUE, offset=0.2, size=tip_fsize[clade_idx], linesize=.3)

p2 = gheatmap(p2, tip_dat, offset=1,
               width=0.4, font.size=3,
               colnames_angle=0,
               colnames_level=reg_lbl,
               colnames_position="top",
               colnames_offset_y=1)

p2 = p2 + scale_fill_manual(name="Regions",
                            breaks=LETTERS[1:5], 
                            values=my_colors,
                            labels=my_labels[[clade_idx]],
                            na.value="transparent")
p2 = p2 + theme(
  legend.position = "right",
  legend.justification = "top"
)

p2


pdf(file=plot_fn[[clade_idx]], height=12, width=7)
print(p2)
dev.off()


