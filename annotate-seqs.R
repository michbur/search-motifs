library(dplyr)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

all_prots <- "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/raw_seq.fasta"

# SignalP 4.1 --------------------------------------

signalP_command <- paste0(files_path, "signalp-4.1/signalp -t gram- -f short ",
                          all_prots, " > signalP.txt")
system(signalP_command)
res_signalP <- read.table("signalP.txt", skip = 1, stringsAsFactors = FALSE) %>%
  select(name = V1,
         signalP_prob = V9,
         signalP_cs = V5) %>%
  mutate(signalP_cs = signalP_cs - 1)


# SignalP 5.0 ---------------------------------------------
signalP5_command <- paste0(files_path, "signalP5/bin/signalp-5 -organism Gram-negative -params_dir '",
                           files_path, "signalP5/bin/model_params' -fasta ", all_prots)

system(signalP5_command)
res_signalP5 <- read.table("output_prediction_summary.txt", stringsAsFactors = FALSE, sep = "\t") %>%
  select(name = V1, signalP5_prob = V3, signalP5_cs = V7) %>%
  mutate(signalP5_prob = as.numeric(signalP5_prob),
         signalP5_cs = as.character(signalP5_cs),
         signalP5_cs = sapply(strsplit(signalP5_cs, ".", fixed = TRUE), dplyr::first),
         signalP5_cs = sapply(strsplit(signalP5_cs, " "), last),
         signalP5_cs = sapply(strsplit(signalP5_cs, "-"), dplyr::first),
         signalP5_cs = as.numeric(signalP5_cs))

# DeepSig ----------------------------


deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ",
                          all_prots, " gramn deepSig.txt")

system(deepSig_command)
res_deepSig <- read.table("deepSig.txt", stringsAsFactors = FALSE) %>%
  mutate(deepSig_prob = ifelse(V2 == "SignalPeptide", V3, 1 - V3),
         deepSig_cs = as.numeric(ifelse(V4 == "-", NA, V4))) %>%
  select(name = V1,
         deepSig_prob,
         deepSig_cs)

# TMHMM ----------------------------

TMHMM_command <- paste0(files_path, "tmhmm-2.0c/bin/tmhmm ", all_prots, "> TMHMM.out")
system(TMHMM_command)
res_TMHMM <- read.table("TMHMM.out", stringsAsFactors = FALSE) %>% 
  select(name = V1, n_hel = V5, topology = V6) %>% 
  mutate(n_hel = as.numeric(sub("PredHel=", "", n_hel)),
         topology = sub("Topology=", "", topology)) 

#file.remove(c("deepSig.txt", "signalP.txt", "output_prediction_summary.txt", "TMHMM.out"))

inner_join(res_signalP, res_deepSig) %>%
  inner_join(res_signalP5) %>%
  inner_join(res_TMHMM) %>%
  dplyr::rename(protein_id = name) %>%
  data.table::fwrite(file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract_annotation.csv", row.names = FALSE,
                     append = FALSE)
