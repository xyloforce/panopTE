library(ggplot2)
library(stringr)
library(zoo)

args = c("data_hg/gc_tes_full_fam_0_0.tsv", "data_pig/gc_tes_full_fam_0_0.tsv", "data_mus/gc_tes_full_fam_0_0.tsv", "data_danio/gc_tes_full_fam_0_0.tsv", "results_hg/savestate_nieb_on_te_full_fam.csv", "results_mus/savestate_nieb_on_te_full_fam.csv", "results_pig/savestate_nieb_on_te_full_fam.csv", "results_danio/savestate_nieb_on_te_full_fam.csv")

source("~/setThemePoster.R")

fix_gc = function(df) {
  df = aggregate(df[, c("count_gc", "total")], by = list(df$te_name, df$pos), FUN = sum)
  df$count_gc = rollapply(df$count_gc,
                          width = 10,
                          FUN = sum,
                          na.rm = TRUE,
                          fill = 0)
  df$total = rollapply(df$total,
                       width = 10,
                       FUN = sum,
                       na.rm = TRUE,
                       fill = 0)
  df$value = df$count_gc / df$total
  colnames(df) = c("te_name", "pos", "count_gc", "total", "value")
  return(df)
}


col_file = c("te_name", "pos", "strand", "value", "count_gc", "total")
tmp_gc = lapply(args[1:4], function(filename) {
  print(filename)
  tmp = read.delim(filename, header = F, col.names = col_file)
  head(tmp)
  tmp = fix_gc(tmp)
  head(tmp)
  species = str_match(filename, "data_(.*?)/")[, 2]
  tmp$species = species
  return(tmp)
})

data = do.call(rbind, tmp_gc)
data = data[, c("species", "te_name", "pos", "value")]
data$type = "GC"

tmp_niebs = lapply(args[5:8], function(filename) {
  print(filename)
  tmp = read.csv2(filename)
  species = str_match(filename, "results_(.*?)/")[, 2]
  tmp$species = species
  tmp
})
tmp_niebs = do.call(rbind, tmp_niebs)
tmp_niebs = tmp_niebs[, c("species", "te_name", "pos", "rel")]
colnames(tmp_niebs) = c("species", "te_name", "pos", "value")
tmp_niebs$type = "NIEBs"

data = rbind(data, tmp_niebs)
data$species = factor(data$species, levels = c("hg", "pig", "mus", "danio"))
# data = rbind(data1, data2, data3, data4)

plot = ggplot(data = data[data$te_name == "L1" &
                          data$pos %in% -6000:50, ],
              aes(x = pos, y = value)) +
  facet_grid(cols = vars(species),
             rows = vars(type),
             scales = "free_y",
             labeller = labeller(species = c("hg" = "H. sapiens",
                                             "mus" = "M. musculus",
                                             "pig" = "S. scrofa",
                                             "danio" = "D. rerio"))) +
  geom_point() +
  theme_poster +
  xlab("Position from the TE border") +
  ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_L1.png", width = 16, height = 5)

plot = ggplot(data = data[(data$te_name == "Alu" |
                           data$te_name == "tRNA") &
                          data$pos %in% -300:50, ],
              aes(x = pos, y = value)) +
  facet_grid(cols = vars(species),
             rows = vars(type),
             scales = "free_y",
             labeller = labeller(species = c("hg" = "H. sapiens",
                                             "mus" = "M. musculus",
                                             "pig" = "S. scrofa",
                                             "danio" = "D. rerio"))) +
  geom_point() +
  theme_poster +
  xlab("Position from the TE border") +
  ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_Alu.png", width = 16, height = 5)

plot = ggplot(data = data[data$te_name == "hAT-Charlie" &
                          data$pos %in% - 200:50, ],
              aes(x = pos, y = value)) +
  facet_grid(cols = vars(species),
             rows = vars(type),
             scales = "free_y",
             labeller = labeller(species = c("hg" = "H. sapiens",
                                             "mus" = "M. musculus",
                                             "pig" = "S. scrofa",
                                             "danio" = "D. rerio"))) +
  geom_point() +
  theme_poster +
  xlab("Position from the TE border") +
  ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_DNA.png", width = 16, height = 5)

plot = ggplot(data = data[data$te_name == "ERV1" &
                          data$pos %in% - 700:50, ],
              aes(x = pos, y = value)) +
  facet_grid(cols = vars(species),
             rows = vars(type),
             scales = "free_y",
             labeller = labeller(species = c("hg" = "H. sapiens",
                                             "mus" = "M. musculus",
                                             "pig" = "S. scrofa",
                                             "danio" = "D. rerio"))) +
  geom_point() +
  theme_poster +
  xlab("Position from the TE border") +
  ylab("Relative coverage") + labs(color = "Type")
ggsave("plot_gc_ERV1.png", width = 16, height = 5)
