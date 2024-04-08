args = commandArgs(trailingOnly = TRUE)
data = read.delim(args[1],
                  col.names = c("chr", "start", "end",
                                "name", "score", "strand"),
                  header = FALSE)

data$name = paste0(data$chr, ":", data$start, "-",
                 data$end, "(", data$strand, ")")

write.table(data,
            file = args[2],
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t")