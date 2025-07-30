library(ggplot2)
library(zoo)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

print("loading datasets")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("fam", "pos", "strand",
                                "rate", "count", "total"))
data = aggregate(data[, c("count", "total")],
                 by = list(data$fam, data$pos),
                 FUN = sum)
colnames(data) = c("fam", "pos", "count", "total")

window_size = 10
if (length(args) > 4) {
  window_size = as.numeric(args[5])
}

families_list = unique(data$fam)
# combinations = paste(families_list, "+")
# combinations = c(combinations, paste(families_list, "-"))

print("smoothing")
head(data)
# total_combinations = paste(data$fam, data$strand)
final_data = lapply(unique(data$fam), function(fam) {
  subset = data[data$fam == fam, ]
  min_sub = min(subset$pos)
  max_sub = max(subset$pos)
  tmp_df = data.frame(fam = fam,
                      pos = min_sub:max_sub,
                      count = 0,
                      total = 0)
  tmp_df[match(subset$pos, tmp_df$pos), "count"] = subset$count
  tmp_df[match(subset$pos, tmp_df$pos), "total"] = subset$total
  tmp_df$sumC = rollapply(tmp_df$count,
                          window_size,
                          sum,
                          na.rm = TRUE,
                          fill = NA)
  tmp_df$sumT = rollapply(tmp_df$total,
                          window_size,
                          sum,
                          na.rm = TRUE,
                          fill = NA)
  return(tmp_df)
})
final_data = do.call(rbind, final_data)
head(final_data)
final_data$gc = final_data$sumC / final_data$sumT

print("plotting")
write.csv2(final_data,
           file = args[4],
           quote = FALSE)

# colnames(final_data) = c("fruit", "legume", "tondeuse", "hamster", "fabien", "Images", "Musique", "oscar", "tranium")
# summary(final_data)
# plot = ggplot(data = final_data[final_data$pos %in% -50:300, ],
#               aes(x = pos, y = gc, color = strand)) +
#   geom_line() +
#   facet_wrap(~ fam) +
#   geom_vline(xintercept = c(0, 133, 266)) +
#   geom_hline(yintercept = as.numeric(args[2]), color = "red", linewidth = 1) +
#   theme_poster
# ggsave(args[3], height = 20, width = 30)

if (!dir.exists(args[3])) {
  dir.create(args[3])
} else {
  print("directory exists, skipping")
}


for (fam in unique(data$fam)) {
  print(fam)
  plot = ggplot(data = final_data[final_data$pos %in% -50:300 &
                                  final_data$fam == fam, ],
              aes(x = pos, y = gc)) +
  geom_line(linewidth = 2) +
  geom_vline(xintercept = c(0, 133, 266)) +
  geom_hline(yintercept = as.numeric(args[2]), color = "red", linewidth = 1) +
  theme_poster +
  xlab("Position from the NIEB border") +
  ylab("GC content")
  fam = str_replace(fam, "/", "_")
  ggsave(paste0(args[3], "/", fam, ".png"), height = 7, width = 7)
}