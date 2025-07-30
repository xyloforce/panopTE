library(ggplot2)
library(stringr)
library(zoo)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

mean10pb = function(x, n = 10) {
  library(zoo)
  rollapply(x, n, mean, na.rm = TRUE, fill = NA)
}

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

tmp_df = lapply(unique(data$fam), function(fam) {
  subset = data[data$fam == fam, ]
  final_tmp = lapply(unique(subset$type), function(type) {
    subsubset = subset[subset$type == type, ]
    new_df = data.frame("fam" = fam,
                        "type" = type,
                        "pos" = min(subsubset$pos):max(subsubset$pos),
                        "count" = 0,
                        "count_bases" = 0)
    new_df[match(subsubset$pos, new_df$pos), "count"] = subsubset$count
    new_df[match(subsubset$pos, new_df$pos), "count_bases"] = subsubset$count_bases
    new_df$smooth_count = mean10pb(new_df$count, 10)
    new_df$smooth_bases = mean10pb(new_df$count_bases, 10)
    return(new_df)
  })
  final_tmp = do.call(rbind, final_tmp)
  return(final_tmp)
})
data = do.call(rbind, tmp_df)
data$rel = data$smooth_count / data$smooth_bases

write.csv2(data,
           file = args[3],
           quote = FALSE)

if (!dir.exists(args[2])) {
  dir.create(args[2])
} else {
  print("directory exists, skipping")
} #5a008b

for (fam in unique(data$fam)) {
  print(fam)
  plot = ggplot(data = data[data$pos %in% -50:300 &
                            data$fam == fam, ],
                aes(x = pos, y = rel, color = type)) +
              geom_line(linewidth = 1) +
              geom_vline(xintercept = c(0, 133, 266)) +
              theme_poster +
              xlab("Position from the NIEB border") +
              ylab("Substitution rate")
  fam = str_replace(fam, "/", "_")
  ggsave(paste0(args[2], "/", fam, ".png"), height = 10, width = 10)
}