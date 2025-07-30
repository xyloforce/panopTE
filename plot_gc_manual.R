library(ggplot2)

source("~/setThemePoster.R")
count_colnames = c("type", "pos", "strand", "gc", "countGC", "total")

args = c("")

data = read.delim("data_hg/gc_tes_full_fam_0_0.tsv",
                  header = F,
                  col.names = count_colnames)

ggplot(data = data[data$pos %in% -300:100 & data$type == "Alu", ],
       aes(x = pos, y = gc, color = strand)) +
  geom_line() +
  theme_poster
ggsave("gc_Alu_human.png", width = 16, height = 5)

ggplot(data = data[data$pos %in% -6000:1000 & data$type == "L1", ],
       aes(x = pos, y = gc, color = strand)) +
  geom_line() +
  theme_poster
ggsave("gc_L1_human.png")

ggplot(data = data[data$pos %in% -500:100 & data$type == "L2a", ],
       aes(x = pos, y = gc, color = strand)) +
  geom_line() +
  theme_poster
ggsave("gc_L2a_human.png")