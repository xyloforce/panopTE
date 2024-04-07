args = commandArgs(trailingOnly = TRUE)
bed_small_colnames = c("chr", "start", "end")
bed_big_colnames = c(bed_small_colnames, "name", "score", "strand")

polyA = read.delim(args[1], header = FALSE, col.names = bed_big_colnames)
tes = read.delim(args[2], header = FALSE, col.names = bed_big_colnames)

tes$alt_id = paste0(tes$chr, ":", tes$start, "-", tes$end, "(", tes$strand, ")")
polyA$name = tes[match(polyA$name, tes$alt_id), "fam"]
write.table(polyA, file = args[4],
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE)