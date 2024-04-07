library(zoo)
library(stringr)

print("conforming muts...")
args = commandArgs(trailingOnly = TRUE)

bases = read.delim(args[1], header = FALSE,
                   col.names = c("pos", "base", "strand", "count"))

bases = aggregate(bases$count, by = list(bases$pos, bases$base), FUN = sum)
colnames(bases) = c("pos", "base", "count")

muts = read.delim(args[2], header = FALSE,
                  col.names = c("pos", "muts", "strand", "count"))
muts = aggregate(muts$count, by = list(muts$pos, muts$muts), FUN = sum)
colnames(muts) = c("pos", "muts", "count")

muts[, c("base", "bdest")] = str_split_fixed(muts$muts,
                                                    pattern = "_", n = 2)

for (base in unique(muts$base)) {
    subset = bases[bases$base == base, ]
    muts[muts$base == base, "total"] =
        subset[match(muts[muts$base == base, "pos"], subset$pos), "count"]
}

muts = muts[, c("pos", "base", "bdest", "count", "total")]
for (base in unique(muts$base)) {
    for (bdest in unique(muts$bdest)) {
        indexes = which(muts$base == base &
                      muts$bdest == bdest)
        if (length(indexes) > 10) {
            muts[indexes, "sumC"] = rollapply(muts[indexes, "count"],
                                          10, sum, na.rm = TRUE,
                                          fill = "extend")
            muts[indexes, "sumT"] = rollapply(muts[indexes, "total"],
                                          10, sum, na.rm = TRUE,
                                          fill = "extend")
        } else {
            muts[indexes, "sumC"] = sum(muts[indexes, "count"], na.rm = TRUE)
            muts[indexes, "sumT"] = sum(muts[indexes, "total"], na.rm = TRUE)
        }
    }
}

muts$prob = muts$sumC / muts$sumT

muts = muts[, c("pos", "base", "bdest", "prob")]
write.table(muts,
            file = args[3],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)