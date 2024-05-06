library(ggplot2)
library(zoo)
args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")
norm = 10

data = read.delim(args[1], header = FALSE,
                  col.names = c("fam", "pos", "mut",
                                "count", "count_bases"))

data$type = "constant"
data[grep("[AT]_[GC]", data$mut), "type"] = "W → S"
data[grep("[GC]_[AT]", data$mut), "type"] = "S → W"
head(data)
data = aggregate(data[, c("count", "count_bases")],
                 by = list(data$fam, data$pos, data$type),
                 FUN = sum, na.rm = TRUE)
colnames(data) = c("fam", "pos", "type", "count", "count_bases")
data = data[!(is.na(data$count) | is.na(data$count_bases)), ]
data$rel = data$count / data$count_bases

write.csv2(data,
           file = args[3],
           quote = FALSE)
plot = ggplot(data = data[data$pos %in% -50:300, ],
              aes(x = pos, y = rel, color = type)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~fam, scales = "free") +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[2], height = 20, width = 30)

