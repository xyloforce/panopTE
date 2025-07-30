library(zoo)
library(stringr)

print("conforming muts...")
args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1], header = FALSE,
                   col.names = c("id", "pos", "mut", "count", "total"))

final_df = lapply(unique(data$id), FUN = function(id) {
    sub_df = lapply(unique(data$mut), function(mut) {
        min_data = min(data$pos)
        max_data = max(data$pos)
        df = data.frame(pos = min_data:max_data, id = id, mut = mut)
        indexes = which(data$mut == mut & data$id == id)
        df$count = data[indexes, ][match(df$pos, data[indexes, "pos"]), "count"]
        df$total = data[indexes, ][match(df$pos, data[indexes, "pos"]), "total"]
        df$count = rollapply(df$count, 10, sum, na.rm = TRUE, fill = "extend")
        df$total = rollapply(df$total, 10, sum, na.rm = TRUE, fill = "extend")
        df$rate = df$count / df$total
        return(df[, c("pos", "id", "mut", "rate")])
    })
    sub_df = do.call(rbind, sub_df)
})
final_df = do.call(rbind, final_df)

final_df[, c("base", "bdest")] =
    str_split_fixed(final_df$mut, pattern = "_", n = 2)

final_df = final_df[, c("id", "pos", "base", "bdest", "rate")]
write.table(final_df,
            file = args[2],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)