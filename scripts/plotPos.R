library(ggplot2)
library(stringr)
library(doParallel)
library(foreach)
library(zoo)

mean10pb = function(x, n = 10) {
  library(zoo)
  rollapply(x, n, mean, na.rm = TRUE, fill = NA)
}

args = commandArgs(trailingOnly = TRUE)

registerDoParallel(cores = as.numeric(args[6]))
source("~/setThemePoster.R")

print("loading data")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "strand",
                                "strand_2", "pos", "val"))
ref = read.delim(args[2], header = FALSE, col.names = c("pos", "val"))

data$strand_m = paste0(data$strand, data$strand_2)
data[data$strand_m %in% c("++", "--"), "strand_m"] = "+"
data[data$strand_m != "+", "strand_m"] = "-"

print("aggregating")
agg_data = aggregate(data$val,
                     by = list(data$id, data$pos, data$strand_m),
                     FUN = sum)
colnames(agg_data) = c("type", "pos", "strand", "value")

total_bases = aggregate(data$val,
                        by = list(data$id, data$strand_m),
                        FUN = sum)
colnames(total_bases) = c("type", "strand", "total")

selected_tes =
  intersect(total_bases[total_bases$strand == "+" &
                        total_bases$total > 20, "type"],
            total_bases[total_bases$strand == "-" &
                        total_bases$total > 20, "type"])
# make sure we have a start and end point for foreach loop
head(selected_tes)

print("cleaning the signal")
result = foreach(name = selected_tes) %dopar% {
  plus = agg_data[agg_data$type == name &
                  agg_data$strand == "+", ]
  tmp = data.frame(
    type = name,
    strand = "+",
    pos = seq(min(plus$pos, na.rm = T),
              max(plus$pos, na.rm = T),
              by = 1),
    value = 0
  )
  tmp[match(plus$pos, tmp$pos), "value"] = plus$value
  plus = tmp
  plus$smooth = mean10pb(plus$value)

  minus = agg_data[agg_data$type == name &
                   agg_data$strand == "-", ]
  tmp = data.frame(
    type = name,
    strand = "-",
    pos = seq(min(minus$pos, na.rm = T),
              max(minus$pos, na.rm = T),
              by = 1),
    value = 0
  )
  tmp[match(minus$pos, tmp$pos), "value"] = minus$value
  minus = tmp
  minus$smooth = mean10pb(minus$value)
  return(rbind(plus, minus))
}

result = do.call(rbind, result)
ids = paste(result$type, result$strand)
ids_ref = paste(total_bases$type, total_bases$strand)

result$rel = result$smooth /
  ref[match(result$pos, ref$pos), "val"]

result$rel = result$rel /
  (total_bases[match(ids, ids_ref), "total"] / sum(ref$val))

write.csv2(result,
           file = args[4],
           quote = FALSE)

result$is_in = result$pos %in% -100:600 # check if is in the plot
counts_per_type = aggregate(result$value, by = list(result$type, result$is_in), FUN = sum)
colnames(counts_per_type) = c("type", "in_plot", "value")
ref$is_in = ref$pos %in% -100:600
ref_per_type = aggregate(ref$val, by = list(ref$is_in), FUN = sum)
colnames(ref_per_type) = c("in_plot", "value")
counts_per_type$ref = ref_per_type[match(counts_per_type$in_plot,
                                         ref_per_type$in_plot),
                                   "value"]
write.csv2(counts_per_type, file = args[5])

if (!dir.exists(args[3])) {
  dir.create(args[3])
} else {
  print("directory exists, skipping")
} #5a008b

print("plotting")
for (type in unique(result$type)) {
  if (substr(type, 1, 1) != "(") {
    print(type)
    plot = ggplot(data = result[result$pos %in% -100:600 &
                                result$type == type, ],
                  aes(x = pos, y = rel, color = strand)) +
      geom_line(linewidth = 1) +
      # geom_hline(yintercept = 1, color = "darkred") +
      labs(x = NULL, y = NULL) +
      theme_poster +
      theme(legend.position = "none") +
      geom_vline(xintercept = c(0, 133, 266),
                 size = 1,
                 color = "darkred")
    type = str_replace(type, "/", "_")
    ggsave(paste0(args[3], "/", type, ".png"), height = 7, width = 7)
  } # remove all single repeats that are SOOOO NUMEROUS
}
# ggtitle(paste(type, sum(result[result$pos %in% -100:600 & result$type == type, "value"]), sep = "/"))