bed_small_colnames = c("chr", "start", "end")
bed_big_colnames = c(bed_small_colnames, "name", "score", "strand")

data = read.delim("results/niebsXtes_shuffled.tsv", header = FALSE,
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
data$id_niebs = make_id(data)
data$id_tes = make_id(data, alt_cols)
data$te_len = data$end.1 - data$start.1

family = read.delim("data/dfam_families.tsv", header = FALSE,
                    col.names = c("name", "family", "subfamily",
                                  "len", "count"))

te_annot = read.delim("data/tes_normal.bed", header = FALSE,
                      col.names = bed_big_colnames)
total_annot = table(te_annot$name)
total_annot = data.frame("name" = names(total_annot),
                         "count" = as.vector(total_annot))
family$count = total_annot[match(family$name, total_annot$name), "count"]

data$family = family[match(data$name, family$name), "family"]
pos = match(data$name, family$name)
min_v = family[pos, "len"] - 0.1 * family[pos, "len"]
max_v = family[pos, "len"] + 0.1 * family[pos, "len"]
data$complete = data$te_len < max_v & data$te_len > min_v
fam_by_status = aggregate(data$family,
                          by = list(data$family, data$complete),
                          FUN = length)
colnames(fam_by_status) = c("fam", "complete", "count")
fam_total = aggregate(family$count,
                      by = list(family$family),
                      FUN = sum)
colnames(fam_total) = c("family", "count")
fam_by_status$relative = fam_by_status$count /
        fam_total[match(fam_by_status$fam, fam_total$family), "count"]

ggplot(data = fam_by_status, aes(x = fam, y = relative, fill = complete)) +
                geom_col(position = "dodge")
# head(data[is.na(data$family),])
# kept_names = names(sort(table(data$name), decreasing = TRUE)[1:100])

# ggplot(data = data[data$name %in% kept_names,], aes(x = te_len)) + geom_histogram() + geom_vline(data = family[family$name %in% kept_names, ], aes(xintercept = len)) + facet_wrap(~name) + xlim(NA, max(family[family$name %in% kept_names, "len"]))