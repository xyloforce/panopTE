library(ggplot2)
library(stringr)

source("~/setThemePoster.R")
# args = c("equilibrium_all_hg.tsv", "equilibrium_all_pig.tsv", "equilibrium_all_mus.tsv", "equilibrium_all_danio.tsv", "results_hg/savestate_nuc_on_te_full_fam.csv", "results_mus/savestate_nuc_on_te_full_fam.csv", "results_pig/savestate_nuc_on_te_full_fam.csv", "results_danio/savestate_nuc_on_te_full_fam.csv")
args = c("equilibrium_all_hg.tsv", "equilibrium_all_pig.tsv", "equilibrium_all_mus.tsv", "equilibrium_all_danio.tsv", "results_hg/savestate_nieb_on_te_full_fam.csv", "results_mus/savestate_nieb_on_te_full_fam.csv", "results_pig/savestate_nieb_on_te_full_fam.csv", "results_danio/savestate_nieb_on_te_full_fam.csv")

tmp_equilibrium = lapply(args[1:4], function(filename) {
  print(filename)
  tmp = read.delim(filename, col.names = c("id", "pos", "mean"))
  species = str_match(filename, "equilibrium_[a-z]+_(.*?)\\.")[, 2]
  tmp$species = species
  tmp
})
tmp_equilibrium = do.call(rbind, tmp_equilibrium)

tmp_niebs = lapply(args[5:8], function(filename) {
  print(filename)
  tmp = read.csv2(filename)
  species = str_match(filename, "results_(.*?)/")[, 2]
  tmp$species = species
  tmp
})
tmp_niebs = do.call(rbind, tmp_niebs)
tmp_niebs = tmp_niebs[, c("species", "te_name", "pos", "rel")]
colnames(tmp_niebs) = c("species", "id", "pos", "mean")

tmp_equilibrium$pos = as.numeric(tmp_equilibrium$pos)
indexes = tmp_equilibrium$id == "Alu" & tmp_equilibrium$species == "hg"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 300
indexes = tmp_equilibrium$id == "tRNA" & tmp_equilibrium$species == "pig"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 250
indexes = tmp_equilibrium$id == "Alu" & tmp_equilibrium$species == "mus"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 150
indexes = tmp_equilibrium$id == "tRNA" & tmp_equilibrium$species == "danio"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 200

indexes = tmp_equilibrium$id == "L1" & tmp_equilibrium$species == "hg"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 6000
indexes = tmp_equilibrium$id == "L1" & tmp_equilibrium$species == "pig"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 150
indexes = tmp_equilibrium$id == "L1" & tmp_equilibrium$species == "mus"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 1000

tmp_equilibrium[tmp_equilibrium$id == "ERV1", "pos"] = tmp_equilibrium[tmp_equilibrium$id == "ERV1", "pos"] - 300

indexes = tmp_equilibrium$id == "hAT-Charlie" & tmp_equilibrium$species %in% c("pig", "hg", "mus")
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 200
indexes = tmp_equilibrium$id == "hAT-Charlie" & tmp_equilibrium$species == "danio"
tmp_equilibrium[indexes, "pos"] = tmp_equilibrium[indexes, "pos"] - 300

tmp_equilibrium$type = "Equilibrium"
tmp_equilibrium$mean = as.numeric(tmp_equilibrium$mean)
tmp_niebs$mean = as.numeric(tmp_niebs$mean)
factor_equilibrium = mean(tmp_niebs$mean, na.rm = TRUE) / mean(tmp_equilibrium$mean, na.rm = TRUE)
tmp_equilibrium$mean = tmp_equilibrium$mean * factor_equilibrium
tmp_niebs$type = "NIEBs"
# tmp_nucs$type = "nucs"

# data = rbind(tmp_equilibrium[, c("species", "id", "type", "pos", "mean")],
#              tmp_nucs[c("species", "id", "type", "pos", "mean")])
data = rbind(tmp_equilibrium[, c("species", "id", "type", "pos", "mean")],
             tmp_niebs[c("species", "id", "type", "pos", "mean")])
data$mean = as.numeric(data$mean)
data$species = factor(data$species, levels = c("hg", "pig", "mus", "danio"))

plot = ggplot(data[data$id %in% c("Alu", "tRNA") & data$pos %in% -300:100, ], aes(x = pos, y = mean)) + geom_line(size = 1) + facet_grid(cols = vars(species), rows = vars(type), labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio")), scales = "free_y") + theme_poster + xlab("Position from the TE 3' border") + ylab("Mean")
ggsave("equilibrium_nieb_Alu_pretty_loveit_yeah.png", width = 16, height = 5)

plot = ggplot(data[data$id == "L1" & data$pos %in% -6000:100 & data$species == "hg", ], aes(x = pos, y = mean)) + geom_line(size = 1) + facet_grid(cols = vars(species), rows = vars(type), , labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio")), scales = "free_y") + theme_poster + xlab("Position from the TE 3' border") + ylab("Mean")
ggsave("equilibrium_nieb_L1_hg_pretty_loveit_yeah.png", width = 16, height = 5)

plot = ggplot(data[data$id == "L1" & data$pos %in% -1000:100 & data$species %in% c("pig", "mus"), ], aes(x = pos, y = mean)) + geom_line(size = 1) + facet_grid(cols = vars(species), rows = vars(type), , labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio")), scales = "free_y") + theme_poster + xlab("Position from the TE 3' border") + ylab("Mean")
ggsave("equilibrium_nieb_L1_ms_pretty_loveit_yeah.png", width = 16, height = 5)

plot = ggplot(data[data$id == "ERV1" & data$pos %in% -300:100, ], aes(x = pos, y = mean)) + geom_line(size = 1) + facet_grid(cols = vars(species), rows = vars(type), , labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio")), scales = "free_y") + theme_poster + xlab("Position from the TE 3' border") + ylab("Mean")
ggsave("equilibrium_nieb_ERV1_pretty_loveit_yeah.png", width = 16, height = 5)

plot = ggplot(data[data$id == "hAT-Charlie" & data$pos %in% -600:100, ], aes(x = pos, y = mean)) + geom_line(size = 1) + facet_grid(cols = vars(species), rows = vars(type), , labeller = labeller(species = c("hg" = "H. sapiens", "mus" = "M. musculus", "pig" = "S. scrofa", "danio" = "D. rerio")), scales = "free_y") + theme_poster + xlab("Position from the TE 3' border") + ylab("Mean")
ggsave("equilibrium_nieb_hAT-Charlie_pretty_loveit_yeah.png", width = 16, height = 5)
