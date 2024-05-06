library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly=TRUE)

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("id", "pos", "base", "strand", "count"))

data = aggregate(data$count, by = list("id" = data$id,
                                       "pos" = data$pos,
                                       "base" = data$base,
                                       "strand" = data$strand),
                             FUN = sum)
colnames(data) = c("id", "pos", "base", "strand", "count")
# data = read.delim(args[1])
data = data[data$pos < 500, ]
data = data[data$base != "N", ]

## first aggregate to delete strand

data = aggregate(data$count, by = list(data$id, data$pos, data$base), FUN = sum)
colnames(data) = c("id", "pos", "base", "count")

## then aggregate to group GC and non GC

data$GC = "nGC"
data[data$base %in% c("G", "C"), "GC"] = "GC"
data = aggregate(data$count, by = list(data$id, data$pos, data$GC), FUN = sum)
colnames(data) = c("id", "pos", "type", "count")

## then divide GC by sum GC + nGC
# gc = data[data$type == "GC",]
# ngc = data[data$type == "nGC",]
data_gc = reshape(data = data,
                  idvar  = c("id", "pos"),
                  v.names = "count",
                  timevar = "type",
                  direction = "wide")
colnames(data_gc) = c("id", "pos", "count_gc", "count_ngc")
data_gc$total = rowSums(data_gc[, c("count_gc", "count_ngc")], na.rm = TRUE)
data_gc$rate = data_gc$count_gc / data_gc$total

data_gc = data_gc[, c("id", "pos", "rate", "count_gc", "total")]
write.table(data_gc, args[2], row.names = FALSE,
            col.names = FALSE, sep = "\t",
            quote = FALSE)