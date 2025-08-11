library(zoo)
library(stringr)

anyna = function(vector) {
  return(TRUE %in% is.na(vector))
}

print("conforming muts...")
args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1], header = FALSE,
                   col.names = c("id", "pos", "mut", "count", "total"))

final_df = lapply(unique(data$id), FUN = function(id) {
  sub_df = lapply(unique(data$mut), function(mut) {
    indexes = which(data$mut == mut & data$id == id)
    min_data = min(data[indexes, "pos"])
    max_data = max(data[indexes, "pos"])
    if (is.finite(min_data) & is.finite(max_data)) {
      df = data.frame(pos = min_data:max_data,
                      id = id,
                      mut = mut,
                      count = 0,
                      total = 0)
      df[match(data[indexes, "pos"], df$pos), "count"] = data[indexes, "count"]
      df[match(data[indexes, "pos"], df$pos), "total"] = data[indexes, "total"]
      df$count = rollapply(df$count,
                           width = 10,
                           FUN = sum,
                           na.rm = TRUE,
                           fill = 0)
      df$total = rollapply(df$total,
                           width = 10,
                           FUN = sum,
                           na.rm = TRUE,
                           fill = 0)
      df$rate = df$count / df$total
      if (anyna(df$count) | anyna(df$total)) {
        print(head(df[is.na(df$count), ]))
        print(head(df[is.na(df$total), ]))
      }
      df[df$count == 0 & df$total == 0, "rate"] = 0 # fix edge case
      return(df[, c("pos", "id", "mut", "rate")])
    } else {
        print(paste(mut, id))
        return(NULL)
    }
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