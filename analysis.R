library(genbankr)
library(Biostrings)
library(pbapply)
library(dplyr)

source("./SMS/functions.R")


dat <- readGenBank2("./data/UTI89.gbk")
seqs <- cds(dat)


all_cds <- as.vector(mcols(seqs)[["translation"]])

seqs_df <- mcols(seqs)
data.frame(select(data.frame(ranges(seqs)), -width), 
           strand = strand(seqs),
           gene_id = seqs_df[["gene_id"]],
           gene = seqs_df[["gene"]],
           #locus_tag = seqs_df[["locus_tag"]],
           #transcript_id = seqs_df[["transcript_id"]],
           fun = as.list(seqs_df[["function"]]) %>% 
             lapply(paste0, collapse = "; ") %>% 
             unlist,
           protein_id = seqs_df[["protein_id"]],
           product = seqs_df[["product"]],
           note = seqs_df[["note"]])
