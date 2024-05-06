library(ggplot2)
library(zoo)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("fam", "pos", "polyA"))
ref_counts = read.delim(args[2],
                        header = FALSE,
                        col.names = c("pos", "total"))

data$ref = ref_counts[match(data$pos, ref_counts$pos), "total"]
for (fam in unique(data$fam)) {
    indexes = which(data$fam == fam)
    data[indexes, "polyA"] = rollapply(data[indexes, "polyA"],
                                       10, sum, na.rm = TRUE, fill = NA)
    data[indexes, "ref"] = rollapply(data[indexes, "ref"],
                                     10, sum, na.rm = TRUE, fill = NA)
}
data$rel = data$polyA / data$ref
write.csv2(data,
           file = args[4],
           quote = FALSE)
plot = ggplot(data = data[data$pos %in% -50:300, ], aes(x = pos, y = rel)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ fam, scales = "free") +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[3], height = 20, width = 30)