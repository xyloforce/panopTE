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

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "strand",
                                "strand_2", "pos", "val"))
ref = read.delim(args[2], header = FALSE, col.names = c("pos", "val"))
families = read.delim(args[3],
                      header = FALSE,
                      col.names = c("id", "fam", "count", "averagelen"))

data$strand_m = paste0(data$strand, data$strand_2)
data[data$strand_m %in% c("++", "--"), "strand_m"] = "+"
data[data$strand_m != "+", "strand_m"] = "-"
agg_data = aggregate(data$val,
                     by = list(data$id, data$pos, data$strand_m),
                     FUN = sum)
colnames(agg_data) = c("type", "pos", "strand", "value")

total_bases = aggregate(data$val,
                        by = list(data$id, data$strand),
                        FUN = sum)
colnames(total_bases) = c("id", "strand", "total")

total_ref = sum(ref$val)

total_bases$cov = total_bases$total / total_ref

ids = paste(agg_data$type, agg_data$strand)
ids_ref = paste(total_bases$id, total_bases$strand)

tes_count = aggregate(agg_data$type,
                      by = list(agg_data$type,
                                agg_data$strand),
                      FUN = length)
colnames(tes_count) = c("type", "strand", "count")

selected_tes =
  intersect(tes_count[tes_count$strand == "+" &
                      tes_count$count > 10, "type"],
            tes_count[tes_count$strand == "-" &
                      tes_count$count > 10, "type"])

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

result$rel = result$smooth /
  ref[match(result$pos, ref$pos), "val"] /
  total_bases[match(ids, ids_ref), "cov"]

write.csv2(result,
           file = args[5],
           quote = FALSE)

if (!dir.exists(args[4])) {
  dir.create(args[4])
} else {
  print("directory exists, skipping")
} #5a008b

print(head(result[result$pos %in% -50:300, ]))
for (type in unique(result$type)) {
  if (substr(type, 1, 1) != "(") {
    print(type)
    plot = ggplot(data = result[result$pos %in% -50:300 &
                                  result$type == type, ],
                  aes(x = pos, y = smooth, color = strand)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 1, color = "darkred") +
      ylab("Average count") +
      xlab("Position from the NIEB border") +
      theme_poster
    type = str_replace(type, "/", "_")
    ggsave(paste0(args[4], "/", type, ".png"), height = 10, width = 10)
  } # remove all single repeats that are SOOOO NUMEROUS
}
