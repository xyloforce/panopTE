data = read.delim("data/dfam_families.tsv",
                  header = F,
                  col.names = c("name", "fam", "subfam", "size", "count"))
dedup = data[match(unique(data$subfam), data$subfam), c("subfam", "count")]
dedup = dedup[dedup$subfam != "Unknown" & dedup$subfam != "", ]
dedup = dedup[dedup$count > 1000, ]
write.table(dedup$subfam, file = "subfam.txt", col.names = FALSE, row.names = FALSE)