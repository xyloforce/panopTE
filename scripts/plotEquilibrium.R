library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

data = read.delim(args[1])
write.csv2(data,
           file = args[3],
           quote = FALSE)
plot = ggplot(data = data[data$position %in% 0:350, ],
              aes(x = position - 50, y = mean)) +
            geom_line() +
            # geom_col() +
            facet_wrap(~ id) +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
ggsave(args[2], height = 20, width = 30)