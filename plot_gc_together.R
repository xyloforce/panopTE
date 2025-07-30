library(ggplot2)
library(stringr)

args = c("data_hg/gc_tes_full_fam_0_0.tsv", "data_pig/gc_tes_full_fam_0_0.tsv", "data_mus/gc_tes_full_fam_0_0.tsv", "data_danio/gc_tes_full_fam_0_0.tsv")

source("~/setThemePoster.R")

col_file = c("te_name", "pos", "strand", "value", "count_gc", "total")
tmp_gc = lapply(args[1:4], function(filename) {
  print(filename)
  tmp = read.delim(filename, header = F, col.names = col_file)
  tmp = aggregate(tmp$value, by = list(tmp$te_name, tmp$pos), FUN = sum)
  colnames(tmp) = c("te_name", "pos", "value")
  species = str_match(filename, "data_(.*?)/")[, 2]
  tmp$species = species
  return(tmp)
})
data = do.call(rbind, tmp_gc)
data$species = factor(data$species, levels = c("hg", "pig", "mus", "danio"))
# data = rbind(data1, data2, data3, data4)

plot = ggplot(data = data[data$te_name == "L1" & data$pos %in% -6000:50, ], aes(x = pos, y = value)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_L1.png", width = 15, height = 5)

plot = ggplot(data = data[(data$te_name == "Alu" | data$te_name == "tRNA") & data$pos %in% -300:50, ], aes(x = pos, y = value)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_Alu.png", width = 5, height = 5)

plot = ggplot(data = data[(data$te_name == "hAT-Charlie") & data$pos %in% - 200:50, ], aes(x = pos, y = value)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_DNA.png", width = 5, height = 5)

plot = ggplot(data = data[(data$te_name == "ERV1") & data$pos %in% - 700:50, ], aes(x = pos, y = value)) + facet_grid(rows = vars(species), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio"))) + geom_point() + theme_poster + xlab("Position from the TE border") + ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_ERV1.png", width = 10, height = 5)
