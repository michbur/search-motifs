
library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

all_prots <- list.files("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/sequences", full.names = TRUE)

# DeepSig ----------------------------

res_deepSig <- pblapply(all_prots, function(ith_file) {
  deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ",
                            ith_file, " gramn deepSig.txt")
  
  system(deepSig_command)
  read.table("deepSig.txt", stringsAsFactors = FALSE) %>%
    mutate(deepSig_prob = ifelse(V2 == "SignalPeptide", V3, 1 - V3),
           deepSig_cs = as.numeric(ifelse(V4 == "-", NA, V4))) %>%
    select(name = V1,
           deepSig_prob,
           deepSig_cs)
})  %>% 
  bind_rows()

save(res_deepSig, file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/res_deepSig.RData")

# TMHMM ----------------------------

res_TMHMM <-  pblapply(all_prots, function(ith_file) {
  TMHMM_command <- paste0(files_path, "tmhmm-2.0c/bin/tmhmm ", ith_file, "> TMHMM.out")
  system(TMHMM_command)
  read.table("TMHMM.out", stringsAsFactors = FALSE) %>% 
    select(name = V1, n_hel = V5, topology = V6) %>% 
    mutate(n_hel = as.numeric(sub("PredHel=", "", n_hel)),
           topology = sub("Topology=", "", topology)) 
})

save(res_TMHMM, file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/res_TMHMM.RData")

#file.remove(c("deepSig.txt", "signalP.txt", "output_prediction_summary.txt", "TMHMM.out"))

# inner_join(res_signalP, res_deepSig) %>%
#   inner_join(res_signalP5) %>%
#   inner_join(res_TMHMM) %>%
#   dplyr::rename(protein_id = name) %>%
#   data.table::fwrite(file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract_annotation.csv", row.names = FALSE,
#                      append = FALSE)
