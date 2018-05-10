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
# load("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract.RData")
# 
# res_df <- bind_rows(res[sapply(res, function(i) class(i) %in% "data.frame")])
# 
# seqinr::write.fasta(sequences = strsplit(res_df[["translation"]], ""), names = res_df[["protein_id"]],
#                     file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/raw_seq.fasta")
