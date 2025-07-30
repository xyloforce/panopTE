library(ggplot2)
library(stringr)

args = c("results_hg/savestate_nieb_on_te_full.csv", "results_mus/savestate_nieb_on_te_full.csv", "results_pig/savestate_nieb_on_te_full.csv", "results_danio/savestate_nieb_on_te_full.csv", "results_hg/savestate_nuc_on_te_full.csv", "results_mus/savestate_nuc_on_te_full.csv", "results_pig/savestate_nuc_on_te_full.csv", "results_danio/savestate_nuc_on_te_full.csv")

source("~/setThemePoster.R")

tmp_niebs = lapply(args[1:4], function(filename) {
  print(filename)
  tmp = read.csv2(filename)
  species = str_match(filename, "results_(.*?)/")[, 2]
  tmp$species = species
  tmp
})
tmp_niebs = do.call(rbind, tmp_niebs)

tmp_nucs = lapply(args[5:8], function(filename) {
  print(filename)
  tmp = read.csv2(filename)
  species = str_match(filename, "results_(.*?)/")[, 2]
  tmp$species = species
  tmp
})
tmp_nucs = do.call(rbind, tmp_nucs)
tmp_nucs = aggregate(tmp_nucs$rel, by = list(tmp_nucs$species, tmp_nucs$te_name, tmp_nucs$pos), FUN = mean)
colnames(tmp_nucs) = c("species", "te_name", "pos", "rel")

# tmp_niebs$type = "normal"


tmp_niebs$type = "niebs"
tmp_nucs$type = "nucs"
data = rbind(tmp_nucs[, c("species", "te_name", "type", "pos", "rel")],
             tmp_niebs[, c("species", "te_name", "type", "pos", "rel")])

# data = rbind(data1, data2, data3, data4)

plot = ggplot(data = data[data$te_name == "L1MB3" & data$pos %in% - 500:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species)) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_L1MB3.png", width = 16, height = 5)

plot = ggplot(data = data[(data$te_name == "AluJb" | data$te_name == "B1_Mm" | data$te_name == "HE2_DR") & data$pos %in% - 300:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_AluJb_and_associates.png", width = 16, height = 5)

# plot = ggplot(data = data[(data$te_name == "hAT-Charlie") & data$pos %in% - 200:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species)) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type") + scale_y_log10()
# ggsave("plot_nuc_DNA.png", width = 16, height = 5)

plot = ggplot(data = data[(data$te_name == "ERV1") & data$pos %in% - 500:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_ERV1.png", width = 16, height = 5)

plot = ggplot(data = data[data$te_name == "L2a" & data$pos %in% - 500:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_L2a.png", width = 16, height = 5)