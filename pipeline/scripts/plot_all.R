library(ggplot2)
library(zoo)
args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")
norm = 10

# args are gc, polyA, muts, plotgc, plotpolya, plotmuts
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("pos", "countGC", "gc_content",
                                "rate", "count", "total", "fam"))

for (fam in unique(data$fam)) {
    data[data$fam == fam, "sumC"] =
                rollapply(data[data$fam == fam, "count"],
                          10, sum, na.rm = TRUE, fill = NA)
    data[data$fam == fam, "sumT"] =
                rollapply(data[data$fam == fam, "total"],
                          10, sum, na.rm = TRUE, fill = NA)
}
data$gc = data$sumC / data$sumT
summary(data)

plot = ggplot(data = data[data$pos %in% 0:350, ], aes(x = pos - 50, y = gc)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ fam) +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[5], height = 20, width = 30)

data = read.delim(args[2],
                  header = FALSE,
                  col.names = c("fam", "pos", "polyA"))

plot = ggplot(data = data[data$pos %in% -50:300, ], aes(x = pos, y = polyA)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ fam, scales = "free") +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[6], height = 20, width = 30)

agg = 10
data = read.delim(args[3], header = FALSE,
                  col.names = c("fam", "pos", "type",
                                "count", "count_bases", "rate"))
data$sumC = NA
data$sumT = NA
for (fam in unique(data$fam)) {
  for (type in unique(data$type)) {
      indexes = which(data$type == type & data$fam == fam)
      # print(indexes)
      # print(length(indexes))
      # print(nrow(data))
      # print(head(rollapply(data[indexes, "count"],
      #                       10, sum, na.rm = TRUE, fill = NA)))
      # print(head(rollapply(data[indexes, "count_bases"],
      #                       10, sum, na.rm = TRUE, fill = NA)))
      # print(fam)
      # print(type)
      data[indexes, "sumC"] =
                  rollapply(data[indexes, "count"],
                            10, sum, na.rm = TRUE, fill = NA)
      data[indexes, "sumT"] =
                  rollapply(data[indexes, "count_bases"],
                            10, sum, na.rm = TRUE, fill = NA)
  }
}
data$rel = data$sumC / data$sumT
plot = ggplot(data = data[data$pos %in% -50:300, ],
              aes(x = pos, y = rel, color = type)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~fam) +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[7], height = 20, width = 30)

data = read.delim(args[4],
                  header = FALSE,
                  col.names = c("fam", "pos", "total_gc", "gc", "meanGC"))

plot = ggplot(data = data[data$pos %in% 0:350, ], aes(x = pos - 50, y = meanGC)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ fam) +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[8], height = 20, width = 30)