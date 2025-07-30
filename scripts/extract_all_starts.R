options(scipen = 999)
args = commandArgs(trailingOnly = TRUE)

bed_colnames = c("chr", "start", "end",
                 "name", "score", "strand")
data = read.delim(args[1],
                  header = TRUE,
                  col.names = bed_colnames)

# ids = read.table(args[2], header = F)$V1

# ids_data = paste(data$chr,
#                  data$start,
#                  data$end,
#                  sep = "_")
# head(ids_data)
# head(ids)
# data = data[ids_data %in% ids, ]

prev_end = c(NA, data[1:(nrow(data) - 1), "end"])
next_start = c(data[2:nrow(data), "start"], NA)
# get chrs to be sure that we aren't going outside of chromosomal boundaries
prev_chr = c(NA, data[1:(nrow(data) - 1), "chr"])
next_chr = c(data[2:nrow(data), "chr"], NA)

data[, c("prev_end", "next_start")] = 
  cbind(prev_end,
        next_start)
data[, c("prev_chr", "next_chr")] =
  cbind(prev_chr,
        next_chr)

summary(data)

plus = data[data$strand == "+", ]
plus$zero = plus$end
plus$end = round((plus$end + plus$next_start) / 2)
plus = plus[plus$chr == plus$next_chr, ]

# plus$zero = plus$end
# plus$start = plus$zero - 50
# plus$end = plus$end + 50

moins = data[data$strand == "-", ]
moins$zero = moins$start
moins$start = round((moins$start + moins$prev_end) / 2)
moins = moins[moins$chr == moins$prev_chr, ]
# moins$zero = moins$start
# moins$end = moins$end + 50
# moins$start = moins$start - 50

data = rbind(plus, moins)
data = data[!is.na(data$start), ]
data = data[!is.na(data$end), ]

print("writing")

write.table(data[, c(bed_colnames, "zero")],
            file = args[2],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
