library(ggplot2)
library(zoo)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)
registerDoParallel(cores = 20)
source("~/setThemePoster.R")

print("loading datasets")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("fam", "pos", "strand",
                                "rate", "count", "total"))

window_size = 10
if (length(args) > 4) {
  window_size = as.numeric(args[5])
}

families_list = unique(data$fam)
combinations = paste(families_list, "+")
combinations = c(combinations, paste(families_list, "-"))

print("smoothing")
head(data)
total_combinations = paste(data$fam, data$strand)
final_data = foreach(x = combinations) %dopar% {
  indexes = which(total_combinations == x)
  subset = data[indexes, ]
  subset$sumC = rollapply(subset$count,
                          window_size,
                          sum,
                          na.rm = TRUE,
                          fill = NA)
  subset$sumT = rollapply(subset$total,
                          window_size,
                          sum,
                          na.rm = TRUE,
                          fill = NA)
  return(subset)
}
final_data = do.call(rbind, final_data)
head(final_data)
final_data$gc = final_data$sumC / final_data$sumT

print("plotting")
write.csv2(final_data,
           file = args[4],
           quote = FALSE)

# colnames(final_data) = c("fruit", "legume", "tondeuse", "hamster", "fabien", "Images", "Musique", "oscar", "tranium")
summary(final_data)
plot = ggplot(data = final_data[final_data$pos %in% -50:300, ],
              aes(x = pos, y = gc, color = strand)) +
  geom_line() +
  facet_wrap(~ fam) +
  geom_vline(xintercept = c(0, 133, 266)) +
  geom_hline(yintercept = as.numeric(args[2]), color = "red", linewidth = 1) +
  theme_poster
ggsave(args[3], height = 20, width = 30)
