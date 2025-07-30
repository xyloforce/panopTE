library(ggplot2)
library(stringr)

args = c("equilibrium_all_hg.tsv", "equilibrium_all_pig.tsv", "equilibrium_all_mus.tsv", "equilibrium_all_danio.tsv", "results_hg/savestate_nuc_on_te_full_fam.csv", "results_mus/savestate_nuc_on_te_full_fam.csv", "results_pig/savestate_nuc_on_te_full_fam.csv", "results_danio/savestate_nuc_on_te_full_fam.csv")

tmp_equilibrium = lapply(args[1:4], function(filename) {
  print(filename)
  tmp = read.delim(filename, col.names = c("id", "pos", "mean"))
  species = str_match(filename, "equilibrium_(.*?)_")[, 2]
  tmp$species = species
  tmp
})
tmp_equilibrium = do.call(rbind, tmp_equilibrium)

tmp_nucs = lapply(args[5:8], function(filename) {
  print(filename)
  tmp = read.csv2(filename)
  species = str_match(filename, "results_(.*?)/")[, 2]
  tmp$species = species
  tmp
})
tmp_nucs = do.call(rbind, tmp_nucs)
tmp_nucs = tmp_nucs[, c("species", "te_name", "type", "pos", "rel")]
tmp_nucs = tmp_nucs[tmp_nucs$type == "normal", ]
colnames(tmp_nucs) = c("species", "id", "type", "pos", "mean")

tmp_equilibrium[tmp_equilibrium$id == "Alu", "pos"] = tmp_equilibrium[tmp_equilibrium$id == "Alu", "pos"] - 300
tmp_equilibrium[tmp_equilibrium$id == "L1", "pos"] = tmp_equilibrium[tmp_equilibrium$id == "Alu", "pos"] - 6000
tmp_equilibrium[tmp_equilibrium$id == "ERV1", "pos"] = tmp_equilibrium[tmp_equilibrium$id == "Alu", "pos"] - 300
tmp_equilibrium[tmp_equilibrium$id == "hAT-Charlie", "pos"] = tmp_equilibrium[tmp_equilibrium$id == "Alu", "pos"] - 600

tmp_equilibrium$type = "equilibrium"
tmp_nucs$type = "nucs"

data = rbind(tmp_equilibrium[, c("species", "id", "type", "pos", "mean")],
             tmp_nucs[c("species", "id", "type", "pos", "mean")])

ggplot(data[data$id == "Alu", ], aes(x = pos, y = mean)) + geom_point() + facet_wrap(~ species)