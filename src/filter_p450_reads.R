#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

#############
# FUNCTIONS #
#############

WriteP450Reads <- function(x) {
    my_dt <- best_hits[rf == x]
    my_file <- paste0(out_prefix, x, ".txt")
    fwrite(unique(my_dt[, .(target_name)]),
           my_file,
           col.names = FALSE)
}

###########
# GLOBALS #
###########

# tidy read of hmmer tblout
cn <- c("target_name",
        "target_accession",
        "query_name",
        "query_accession",
        "full_E-value",
        "full_score",
        "full_bias",
        "domain_E-value",
        "domain_score",
        "domain_bias",
        "domain_no_exp",
        "domain_no_reg",
        "domain_no_clu",
        "domain_no_ov",
        "domain_no_env",
        "domain_no_dom",
        "domain_no_rep",
        "domain_no_inc",
        "description_of_target")

hmmer_results_files <- snakemake@input[["hmmer_results"]]
out_prefix <- snakemake@params[["out_prefix"]]

# dev
# hmmer_results_files <- list.files("output/hmmer",
#                                   pattern = "rf[[:digit:]].txt",
#                                   full.names = TRUE)
# out_prefix <- 'output/hmmer/P450_read_names_rf'

########
# MAIN #
########

# read the hmmer results
names(hmmer_results_files) <- gsub("rf([[:digit:]]).txt",
                                   "\\1",
                                   basename(hmmer_results_files))
hmmer_results_list <- lapply(hmmer_results_files, function(x)
    fread(paste('grep -v "^#"', x), select = 1:19, col.names = cn))
hmmer_results <- rbindlist(hmmer_results_list, idcol = "rf")

# get the best hmmer result for each target
my_rows <- hmmer_results[, .I[which.max(`full_score`)],
                         by = .(target_name)][, V1]
best_hits <- hmmer_results[my_rows][`full_E-value` <= 0.01]
setorder(best_hits, -full_score)

# make a list for each reading frame
lapply(best_hits[, unique(rf)], WriteP450Reads)

# log
sessionInfo()
