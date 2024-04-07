library(ggplot2)
library(extrafont)
library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)

bed_small_colnames = c("chr", "start", "end")
bed_big_colnames = c(bed_small_colnames, "name", "score", "strand")

print("Loading datasets...")
niebs = read.delim(args[1], header = FALSE, col.names = bed_small_colnames)
polyA = read.delim(args[2], header = FALSE, col.names = bed_big_colnames)
tes = read.delim(args[3], header = FALSE, col.names = bed_big_colnames)
families = read.delim(args[4], header = FALSE,
                      col.names = c("name", "family",
                                    "subfamily", "count", "len"))
inter_tes = read.delim(args[5], header = FALSE,
                      col.names = c(bed_small_colnames,
                                    bed_big_colnames,
                                    "overlap"))
inter_tes_shuffled = read.delim(args[6], header = FALSE,
                         col.names = c(bed_small_colnames,
                                       bed_big_colnames,
                                       "overlap"))
inter_polyA = read.delim(args[7], header = FALSE,
                         col.names = c(bed_small_colnames,
                                       bed_big_colnames,
                                       "overlap"))

make_id = function(dataframe, cols = c("chr", "start", "end")) {
        return(paste(dataframe[, cols[1]],
                     dataframe[, cols[2]],
                     dataframe[, cols[3]],
                     sep = "_"))
}

alt_cols = c("chr.1", "start.1", "end.1")
inter_tes$id_niebs = make_id(inter_tes)
inter_tes$id_tes = make_id(inter_tes, alt_cols)
inter_tes_shuffled$id_niebs = make_id(inter_tes_shuffled)
inter_tes_shuffled$id_tes = make_id(inter_tes_shuffled, alt_cols)

print("plotting total numbers")
data = data.frame("count" = c(nrow(niebs),
                              nrow(tes),
                              length(unique(inter_tes$id_niebs)),
                              length(unique(inter_tes$id_tes)),
                              length(unique(inter_tes_shuffled$id_niebs)),
                              length(unique(inter_tes_shuffled$id_tes))),
                  "type" = c("niebs",
                             "tes",
                             "niebs",
                             "tes",
                             "niebs",
                             "tes"),
                  "source" = c("total",
                               "total",
                               "intersection",
                               "intersection",
                               "random",
                               "random"))
plot = ggplot(data = data, aes(x = type, y = count, fill = source)) +
        geom_col(position = "dodge") +
        theme_poster +
        ylab("Number of tes/niebs")
ggsave("total_vs_intersecting_tes_and_niebs.svg", height = 5, width = 10)
kept_tes = names(table(tes$name)[table(tes$name) > 25])

# add family info to interte
inter_tes$family = families[match(inter_tes$name, families$name), "family"]
# add family info to inter_tes_shuffled
inter_tes_shuffled$family = families[match(inter_tes_shuffled$name, families$name), "family"]

# prepare df for family plot
print("plotting tes by fam")
data = sort(table(inter_tes$family), decreasing = TRUE)
data = data.frame("TE_family" = names(data),
                  "count" = as.vector(data),
                  "type" = "normal")

# add shuffled TEs
tmp = sort(table(inter_tes_shuffled$family), decreasing = TRUE)
data = rbind(data, data.frame("TE_family" = names(tmp),
                              "count" = as.vector(tmp),
                              "type" = "shuffled"))

data$TE_family = factor(x = data$TE_family, levels =
                        unique(data[order(data$count), "TE_family"]))
plot = ggplot(data, aes(x = TE_family, y = count, fill = type)) +
        geom_col(position = "dodge") +
        theme_poster +
        theme(axis.text.x = element_text(angle = 90))
ggsave("count_niebs_fam.svg", width = 10, height = 5)

print("plotting relative numbers")
tes$family = families[match(tes$name, families$name), "family"]

tmp = table(tes$family)
data$relative = as.vector(data$count / tmp[match(data$TE_family, names(tmp))])
plot = ggplot(data, aes(x = TE_family, y = relative)) +
        geom_col() + theme_poster +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("TE family")
ggsave("relative_niebs_fam.svg", width = 10, height = 5)

print("plotting polyA tes")
polyA$name = apply(str_match(polyA$name,
                             "(\\w):(\\d+)-(\\d+)")[, 2:4],
                   MARGIN = 1, paste, collapse = "_")
polyA$name = paste0("chr", polyA$name)
head(polyA)
inter_tes$has_polyA = inter_tes$id_tes %in% polyA$name
print(head(inter_tes))
data = aggregate(inter_tes$id_niebs,
                 by = list(inter_tes$family,
                           inter_tes$has_polyA),
                 FUN = length)
colnames(data) = c("TE_family", "has_polyA", "count")
data$relative = as.vector(data$count / tmp[match(data$TE_family, names(tmp))])
data$TE_family = factor(x = data$TE_family, levels = names(sort(tmp)))
plot = ggplot(data, aes(x = TE_family, y = relative, fill = has_polyA)) +
        geom_col() + theme_poster +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("TE family")
ggsave("haspolyA_niebs_fam.svg", width = 10, height = 5)

print("plotting relative overlap")
total_overlap = aggregate(inter_tes$overlap,
                          by = list(inter_tes$family),
                          FUN = sum)
tes$len = tes$end - tes$start
total_length = aggregate(tes$len, by = list(tes$family), FUN = sum)
relative = total_overlap$x / total_length[match(total_overlap$Group.1,
                                              total_length$Group.1), "x"]
data$rel_overlap = relative[match(total_overlap$Group.1,
                                  data$TE_family)]
plot = ggplot(data, aes(x = TE_family, y = rel_overlap, fill = has_polyA)) +
        geom_col() + theme_poster +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("TE family")
ggsave("overlap_niebs_fam.svg", width = 10, height = 5)