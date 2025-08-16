library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)

print("Loading dataset")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "pos", "base", "strand", "count"))

print("restructuring data")
agg_data = aggregate(data$count,
                     by = list(data$id,
                               data$pos,
                               data$base),
                     FUN = sum)
colnames(agg_data) = c("te_name", "pos", "base", "count")

reshaped_data = reshape(data = agg_data,
                        idvar = c("te_name", "pos"),
                        v.names = "count",
                        timevar = "base",
                        direction = "wide")
reshaped_data[is.na(reshaped_data$count.S), "count.S"] = 0
reshaped_data[is.na(reshaped_data$count.W), "count.W"] = 0

reshaped_data$total = reshaped_data$count.S + reshaped_data$count.W
reshaped_data$rel = reshaped_data$count.S /
  reshaped_data$total

head(reshaped_data)
colnames(reshaped_data) = c("id", "pos", "count_gc", "count_ngc", "total", "rate")
print("Writing dataframe")
data_gc = reshaped_data[, c("id", "pos", "rate", "count_gc", "total")]
write.table(data_gc, args[2], row.names = FALSE,
            col.names = FALSE, sep = "\t",
            quote = FALSE)