library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)

file_paths = str_split(args[1], pattern = ":", simplify = TRUE)[1, ]
file_contents = str_split(args[2], pattern = ":", simplify = TRUE)[1, ]

print(file_paths)

data = lapply(seq_along(file_paths), function(x) {
	tmp = read.csv2(file_paths[x])
	tmp$species = file_contents[x]
	return(tmp)
})

data = do.call(rbind, data)
head(data)

subset = data[data$position %in% -50:300, ]
# subset = subset[!is.na(subset$rel),  ]
count_species = aggregate(subset$species, by = list(subset$id), FUN = function(x) {length(unique(x))})
colnames(count_species) = c("id", "species")
count_species = count_species[count_species$species > 1, ]
print(head(count_species))
subset = subset[subset$id %in% count_species$id, ]

plot = ggplot(data = subset, aes(x = position -50, y = mean, color = species)) + geom_line() + facet_wrap(~ id)
ggsave(args[3], width = 10, height = 10)
