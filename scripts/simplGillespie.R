library(stringr)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)
set.seed(1123)

source("~/second_home/scripts/functions.R")
source("~/setThemePoster.R")

safeFixDf = function(muts, offset, total_len) {
    tryCatch({
        suppressWarnings(fix_df(muts, offset, total_len))
    }, error = function(err){
        return(data.frame())
    })
}

registerDoParallel(cores = as.numeric(args[4]))
number_seqs = 1000
number_gens = 2000
gc_content = 10
offset = 50 # offset applied to be able to see mutation rates inside
len_seqs = 350

print("read files")
muts_nCpG = read.delim(args[1])
muts_nCpG = muts_nCpG[muts_nCpG$pos >= -offset &
                      muts_nCpG$pos <= (len_seqs - offset), ]
muts_CpG = read.delim(args[2])
muts_CpG = muts_CpG[muts_CpG$pos >= -offset + 10 &
                    muts_CpG$pos <= (len_seqs - offset), ]

print("correcting dfs...")
muts_nCpG = safeFixDf(muts_nCpG, -offset, len_seqs - offset)
muts_CpG = safeFixDf(muts_CpG, -offset, len_seqs - offset)

if(nrow(muts_CpG) > 0 & nrow(muts_nCpG) > 0) {
    print("corrected")
    seqs = list()
    for (i in 1:number_seqs) {
        seqs[[i]] = define_seq(len = len_seqs, GC = gc_content)
    }
    print("seq defined")

    gc_as_time = data.frame(time = c(), gc_content = c(), seq = c())
    print("simulating...")
    result = foreach(i = seq_len(length(seqs))) %dopar% {
        mutate_seq_2(seqs[[i]], number_gens, muts_nCpG, muts_CpG, offset)
    } # i is gc : seq : i

    print("simulation finished, treating...")
    for (df in lapply(seq_len(length(result)), FUN = function(x) {
        tmp = result[[x]][[2]]
        tmp$seq = x
        return(tmp)
    })) {
        gc_as_time = rbind(gc_as_time, df)
    }
    seqs = lapply(result, FUN = function(x) {
        return(x[[1]])
    })

    cat("\n")
    print("saving")

    profile = get_gc_by_pos(seqs)
    profile$mean = as.numeric(mean10pb(profile$gc_content))

    write.table(profile,
                file = args[3],
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
} else {
    write.table(data.frame("pos" = 1:len_seqs,
                           "gc" = 0,
                           "total" = 0,
                           "rel" = 0,
                           "mean_rel" = 0),
                file = args[3],
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
}
