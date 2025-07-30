library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)

print("Loading dataset")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "pos", "base", "strand", "count"))

data$type = "nGC"
data[data$base == "S", "type"] = "GC"
data = data[, c("id", "pos", "strand", "type", "count")] # reorder cols

print("Reshaping datafame")
data_gc = reshape(data = data,
                  idvar  = c("id", "pos", "strand"),
                  v.names = "count",
                  timevar = "type",
                  direction = "wide")
# colnames(data_gc) = c("id", "pos", "strand", "count_ngc", "count_gc")
colnames(data_gc)[colnames(data_gc) == "count.GC"] = "count_gc"
colnames(data_gc)[colnames(data_gc) == "count.nGC"] = "count_ngc"

data_gc$total = rowSums(data_gc[, c("count_gc", "count_ngc")], na.rm = TRUE)
data_gc$rate = data_gc$count_gc / data_gc$total

print("Writing dataframe")
data_gc = data_gc[, c("id", "pos", "strand", "rate", "count_gc", "total")]
write.table(data_gc, args[2], row.names = FALSE,
            col.names = FALSE, sep = "\t",
            quote = FALSE)
