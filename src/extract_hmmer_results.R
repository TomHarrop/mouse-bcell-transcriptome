#!/usr/bin/env Rscript

library(data.table)

###########
# GLOBALS #
###########


hmmer_results_files <- list.files("output/hmmer",
                                  pattern = "rf[[:digit:]].txt",
                                  full.names = TRUE)

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
best_hits <- hmmer_results[my_rows]
setorder(best_hits, -full_score)

# all transcripts < cutoff
hmmer_results[full_score > 10, unique(target_name)]





