library(stringr)

args = commandArgs(trailingOnly = TRUE)
bases = read.delim(args[1], header = FALSE,
                   col.names = c("id", "pos", "base", "strand", "count"))

muts = read.delim(args[2], header = FALSE,
                  col.names = c("id", "pos", "muts", "strand", "count"))

muts[, c("ancestral", "current")] = str_split_fixed(muts$muts,
                                                    pattern = "_", n = 2)
# reverseBase = function(bases) {
#     tmp = bases
#     tmp[grep("A", bases)] = "T"
#     tmp[grep("T", bases)] = "A"
#     tmp[grep("G", bases)] = "C"
#     tmp[grep("C", bases)] = "G"
#     return(tmp)
# }

# muts[muts$strand == "-", "ancestral"] = reverseBase(muts[muts$strand == "-", "ancestral"])
# muts[muts$strand == "-", "current"] = reverseBase(muts[muts$strand == "-", "current"])

bases = aggregate(bases$count, by = list(bases$id, bases$pos, bases$base), FUN = sum)
colnames(bases) = c("id", "pos", "base", "count")

id_bases = paste(bases$id, bases$pos, bases$base)
id_muts = paste(muts$id, muts$pos, muts$ancestral)
head(id_bases)
head(id_muts)
head(bases)
muts$total = bases[match(id_muts, id_bases), "count"]

write.table(muts[, c("id", "pos", "muts", "count", "total")],
            file = args[3],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
