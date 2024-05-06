library(ggplot2)
library(zoo)
args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("fam", "pos",
                                "rate", "count", "total"))

window_size = 10
if (length(args) > 3) {
    window_size = as.numeric(args[4])
}

for (fam in unique(data$fam)) {
    data[data$fam == fam, "sumC"] =
                rollapply(data[data$fam == fam, "count"],
                          window_size, sum, na.rm = TRUE, fill = NA)
    data[data$fam == fam, "sumT"] =
                rollapply(data[data$fam == fam, "total"],
                          window_size, sum, na.rm = TRUE, fill = NA)
}
data$gc = data$sumC / data$sumT


write.csv2(data,
           file = args[3],
           quote = FALSE)
plot = ggplot(data = data[data$pos %in% -50:300, ], aes(x = pos, y = gc)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ fam) +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[2], height = 20, width = 30)