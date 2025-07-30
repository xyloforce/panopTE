options(scipen = 999)
args = commandArgs(trailingOnly = TRUE)

bed_colnames = c("chr", "start", "end", "name", "score", "strand")
data = read.delim(args[1], header = FALSE, col.names = bed_colnames)
data$name = paste(data$chr, data$start, data$end, sep = "_")

write.table(data, file = args[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)