library(dplyr)
library(genbankr)
library(Biostrings)
library(pbapply)

source("./SMS/functions.R")

#awk -v n=1 '/^$/{close("out"n);n++;next} {print > "/home/michal/Dropbox/ann-arbor-collab/gut-microbiom-data/out"n}' 42773.ADKP01000001-ADKP01000236.nuc.gbk

res <- pblapply(list.files("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-data", full.names = TRUE), 
                function(file_name)
                  try({
                    dat <- data.frame(file = last(strsplit(x = file_name, split = "/", fixed = TRUE)[[1]]),
                                      suppressMessages(readGenBank2(file_name)), stringsAsFactors = FALSE)
                  }, silent = TRUE)
)

save(res, file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract.RData")

load("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract.RData")

res_df <- bind_rows(res[sapply(res, function(i) class(i) %in% "data.frame")])

data.table::fwrite(x = res_df, 
                   file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract.csv", 
                   row.names = FALSE)

seq_splits <- split(1L:nrow(res_df), ceiling(quantile(1L:nrow(res_df), 1L:20/20)))

lapply(1L:length(seq_splits), function(ith_split) {
  seqinr::write.fasta(sequences = strsplit(res_df[seq_splits[[ith_split]], "translation"], ""), 
                      names = res_df[seq_splits[[ith_split]], "protein_id"],
                      file = paste0("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/sequences/split", 
                                    ith_split, ".fasta"))
})
library(dplyr)

load("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/res_deepSig.RData")
load("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/res_signalP.RData")
load("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/res_TMHMM.RData")

inner_join(res_signalP, res_deepSig) %>%
  inner_join(res_TMHMM) %>%
  dplyr::rename(protein_id = name) %>%
  data.table::fwrite(file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract_annotation.csv", row.names = FALSE,
                     append = FALSE)
