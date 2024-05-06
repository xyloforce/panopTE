library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "strand",
                                "strand_2", "pos", "val"))
ref = read.delim(args[2], header = FALSE, col.names = c("pos", "val"))
data$strand_m = paste0(data$strand, data$strand_2)
data[data$strand_m %in% c("++", "--"), "strand_m"] = "+"
data[data$strand_m != "+", "strand_m"] = "-"
data = aggregate(data$val, by = list(data$id, data$pos, data$strand), FUN = sum)
colnames(data) = c("type", "pos", "strand", "value")
data$rel = data$value / ref[match(data$pos, ref$pos), "val"]
write.csv2(data,
           file = args[4],
           quote = FALSE)
plot = ggplot(data = data[data$pos %in% -50:350, ],
              aes(x = pos, y = rel, color = strand)) +
    geom_line() +
    facet_wrap(~ type, scales = "free")
ggsave(args[3], width = 20, height = 20)