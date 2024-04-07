library(stringr)

args = commandArgs(trailingOnly = TRUE)
bases = read.delim(args[1], header = FALSE,
                   col.names = c("pos", "base", "strand", "count"))

muts = read.delim(args[2], header = FALSE,
                  col.names = c("pos", "muts", "strand", "count"))

muts[, c("ancestral", "current")] = str_split_fixed(muts$muts,
                                                    pattern = "_", n = 2)

bases$type = "W"
bases[bases$base %in% c("C", "G"), "type"] = "S"
bases = aggregate(bases$count, by = list(bases$pos, bases$type), FUN = sum)
colnames(bases) = c("pos", "type", "count")

muts$type = "constant"
muts[muts$muts %in% c("A_G", "A_C", "T_G", "T_A"), "type"] = "W → S"
muts[muts$muts %in% c("C_A", "C_T", "G_A", "G_T"), "type"] = "S → W"
muts = aggregate(muts$count, by = list(muts$pos, muts$type), FUN = sum)
colnames(muts) = c("pos", "type", "count")

# muts = muts[muts$type != "constant", ]
muts$rate = 0
convert = list("S" = "S → W", "W" = "W → S", "SW" = "constant")
for (type in convert) {
    if (type == "constant") {
        subset = aggregate(bases$count, by = list(bases$pos),
                           FUN = sum, na.rm = TRUE)
        colnames(subset) = c("pos", "count")
        muts[muts$type == type,
             "ref_bases"] = subset[match(muts[muts$type == type, "pos"],
                                         subset$pos), "count"]
        muts[muts$type == type,
             "rate"] = muts[muts$type == type, "count"] /
                       subset[match(muts[muts$type == type, "pos"],
                                    subset$pos), "count"]
    } else {
        subset = bases[bases$type == names(convert[match(type, convert)]), ]
        muts[muts$type == type,
             "ref_bases"] = subset[match(muts[muts$type == type, "pos"],
                                         subset$pos), "count"]
        muts[muts$type == type,
             "rate"] = muts[muts$type == type, "count"] /
                       subset[match(muts[muts$type == type, "pos"],
                                    subset$pos), "count"]
    }
}

write.table(muts[, c("pos", "type", "count", "ref_bases", "rate")],
            file = args[3],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
