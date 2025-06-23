args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1], header = TRUE)
data = data[data$milliDel < 100 & data$milliIns < 100, ]
# data = data[data$swScore > 624, ]
# data = data[data$repLeft == 0, ] # nothing left on the annotation = we're at the end
plus = data[data$strand == "+" &
            data$repLeft > -20, ]
total_len = plus$repEnd - plus$repStart
plus = plus[total_len > 100, ]

moins = data[data$strand == "-" &
             data$repStart > -20, ]
total_len = moins$repLeft - moins$repStart
moins = moins[total_len > 100, ]

data = rbind(plus,
             moins)
head(data)
write.table(data,
            row.names = FALSE,
            col.names = TRUE,
            file = args[2],
            quote = FALSE,
            sep = "\t")
