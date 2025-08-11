library(ggplot2)
library(zoo)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

print("loading data")
data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("fam", "strandA", "strandB", "pos", "polyA"))
data$strand = paste0(data$strandA, data$strandB)
data[data$strand %in% c("++", "--"), "strand"] = "+"
data[data$strand %in% c("+-", "-+"), "strand"] = "-"
data = aggregate(data$polyA, by = list(data$fam, data$pos, data$strand), FUN = sum)
colnames(data) = c("fam", "pos", "strand", "polyA")

ref_counts = read.delim(args[2],
                        header = FALSE,
                        col.names = c("pos", "total"))

data$ref = ref_counts[match(data$pos, ref_counts$pos), "total"]

print("smoothing counts")
list_df = lapply(unique(data$fam), function(fam) {
    sublist = lapply(c("-", "+"), function(strand) {
        subset = data[data$fam == fam &
                      data$pos %in% -50:350 &
                      data$strand == strand, ]
        tmp_df = data.frame(fam = fam,
                            pos = -50:350,
                            strand = strand,
                            polyA = 0,
                            ref = 0)
        tmp_df[match(subset$pos, tmp_df$pos), "polyA"] = subset$polyA
        tmp_df[match(subset$pos, tmp_df$pos), "ref"] = subset$ref
        tmp_df$polyA = rollapply(tmp_df$polyA,
                                10, sum,
                                na.rm = TRUE, fill = NA)
        tmp_df$ref = rollapply(tmp_df$ref,
                            10, sum,
                            na.rm = TRUE, fill = NA)
        return(tmp_df)
    })
    sublist = do.call(rbind, sublist)
    return(sublist)
})
data = do.call(rbind, list_df)
data$rel = data$polyA / data$ref
write.csv2(data,
           file = args[4],
           quote = FALSE)

if (!dir.exists(args[3])) {
  dir.create(args[3])
} else {
  print("directory exists, skipping")
} #5a008b

print("plotting")
for (fam in unique(data$fam)) {
  if (substr(fam, 1, 1) != "(") {
    print(fam)
    plot = ggplot(data = data[data$pos %in% -50:350 &
                              data$fam == fam, ],
                  aes(x = pos, y = rel, color = strand)) +
            geom_line() +
            geom_vline(xintercept = c(0, 133, 266)) +
            theme_poster
    ggsave(paste0(args[3], "/", fam, ".png"), height = 7, width = 7)
  } # remove all single repeats that are SOOOO NUMEROUS
}