args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1], header = TRUE)
data = data[data$milliDel < 100 & data$milliIns < 100, ]
data = data[data$milliDel < 100 & data$milliIns < 100, ]
# try
data$len = data$repEnd - data$repStart
data[data$len > 100, ]

plus = data[data$strand == "+" &
            data$repLeft > -20, ]
moins = data[data$strand == "-" &
             data$repStart > -20, ]

data = rbind(plus,
             moins)
head(data)
data = data[!grepl("_", data$genoName), ]
write.table(data,
            row.names = FALSE,
            col.names = TRUE,
            file = args[2],
            quote = FALSE,
            sep = "\t")
