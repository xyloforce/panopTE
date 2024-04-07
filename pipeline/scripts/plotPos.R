library(ggplot2)

data = read.delim("data/pos_full.tsv", header = FALSE)
ref = read.delim("ref_counts.tsv", header = F)
data$strand = paste0(data$V1, data$V2)
data[data$strand %in% c("++", "--"), "strand"] = "+"
data[data$strand != "+", "strand"] = "-"
data = aggregate(data$V4, by = list(data$V3, data$V5, data$strand), FUN = sum)
colnames(data) = c("pos", "type", "strand", "value")
data$rel = data$value / ref[match(data$pos, ref$V1), "V2"]
plot = ggplot(data = data[data$pos %in% -50:350,], aes(x = pos, y = rel, color = strand)) +
    geom_line() +
    facet_wrap(~ type, scales = "free")
ggsave("pos_te_by_fam.png", width = 20, height = 20)