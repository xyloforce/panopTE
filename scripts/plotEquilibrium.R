library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

data = read.delim(args[1])
write.csv2(data,
           file = args[3],
           quote = FALSE)

if (!dir.exists(args[2])) {
  dir.create(args[2])
} else {
  print("directory exists, skipping")
} #5a008b

for (id in unique(data$id)) { # 
  print(id)
  plot = ggplot(data = data[data$position %in% 0:350 &
                            data$id == id, ],
                aes(x = position - 50, y = mean)) +
            geom_line(linewidth = 1) +
            # geom_col() +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster +
            xlab("Position from the NIEB border") +
            ylab("GC content at equilibrium")
  ggsave(paste0(args[2], "/", id, ".png"), height = 10, width = 10)
}