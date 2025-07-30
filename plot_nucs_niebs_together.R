library(ggplot2)
library(stringr)

args = c("results_hg/savestate_nieb_on_te_full_fam.csv", "results_mus/savestate_nieb_on_te_full_fam.csv", "results_pig/savestate_nieb_on_te_full_fam.csv", "results_danio/savestate_nieb_on_te_full_fam.csv", "results_hg/savestate_nuc_on_te_full_fam.csv", "results_mus/savestate_nuc_on_te_full_fam.csv", "results_pig/savestate_nuc_on_te_full_fam.csv", "results_danio/savestate_nuc_on_te_full_fam.csv")

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
data$species = factor(data$species, levels = c("hg", "pig", "mus", "danio"))

plot = ggplot(data = data[data$te_name == "L1" & data$pos %in% - 6000:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_L1.png", width = 15, height = 5)

plot = ggplot(data = data[(data$te_name == "Alu" | data$te_name == "tRNA") & data$pos %in% - 300:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_Alu.png", width = 5, height = 5)

plot = ggplot(data = data[(data$te_name == "hAT-Charlie") & data$pos %in% - 200:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type") + scale_y_log10()
ggsave("plot_nuc_DNA.png", width = 5, height = 5)

plot = ggplot(data = data[(data$te_name == "ERV1") & data$pos %in% - 700:50, ], aes(x = pos, y = rel, color = type)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_nuc_ERV1.png", width = 10, height = 5)
