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

if(length(args) > 4) {
    offset = as.numeric(args[5])
    len_seqs = as.numeric(args[6])
}

print("read files")
whole_muts_nCpG = read.delim(args[1])
whole_muts_CpG = read.delim(args[2])

colnames(whole_muts_CpG) = c("id", "pos", "base", "bdest", "prob")
colnames(whole_muts_nCpG) = c("id", "pos", "base", "bdest", "prob")

print("working on tes...")
final_results = lapply(unique(whole_muts_nCpG$id), FUN = function(id) {
    profile = data.frame()
    print(paste0("id : ", id))
    muts_nCpG = whole_muts_nCpG[whole_muts_nCpG$id == id, ]
    muts_CpG = whole_muts_CpG[whole_muts_CpG$id == id, ]
    muts_nCpG = safeFixDf(muts_nCpG, -offset, len_seqs - offset)
    muts_CpG = safeFixDf(muts_CpG, -offset, len_seqs - offset)

    if (nrow(muts_CpG) > 0 && nrow(muts_nCpG) > 0) {
        seqs = list()
        for (i in 1:number_seqs) {
            seqs[[i]] = define_seq(len = len_seqs, GC = gc_content)
        }

        gc_as_time = data.frame(time = c(), gc_content = c(), seq = c())
        print("simulating...")
        result = foreach(i = seq_len(length(seqs))) %dopar% {
            mutate_seq_2(seqs[[i]], number_gens, muts_nCpG, muts_CpG, offset)
        } # i is gc : seq : i

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

        profile = get_gc_by_pos(seqs)
        profile$mean = as.numeric(mean10pb(profile$gc_content))
        profile$id = id
        profile = profile[, c("id", "position", "mean")]
    }
    return(profile)
})

print("saving results...")
final_results = do.call(rbind, final_results)

write.table(final_results,
            file = args[3],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
