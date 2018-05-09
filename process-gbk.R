library(dplyr)
library(genbankr)
library(Biostrings)
library(pbapply)

#awk -v n=1 '/^$/{close("out"n);n++;next} {print > "/home/michal/Dropbox/ann-arbor-collab/gut-microbiom-data/out"n}' 42773.ADKP01000001-ADKP01000236.nuc.gbk

res <- pblapply(list.files("/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-data/", full.names = TRUE), function(file_name) 
  try({
    
    dat <- suppressMessages(readGenBank(file_name))
    
    seqs <- cds(dat)
    
    all_cds <- as.vector(mcols(seqs)[["translation"]])
    
    seqs_df <- mcols(seqs)
    gene_fun <- as.list(seqs_df[["function"]]) %>% 
      lapply(paste0, collapse = "; ") %>% 
      unlist

    df_list <- list(source = definition(dat),
                    file = last(strsplit(file_name, "/")[[1]]),
                    select(data.frame(ranges(seqs)), -width), 
                    strand = strand(seqs),
                    gene_id = seqs_df[["gene_id"]],
                    gene = seqs_df[["gene"]],
                    locus_tag = seqs_df[["locus_tag"]],
                    transcript_id = seqs_df[["transcript_id"]],
                    fun = gene_fun,
                    protein_id = seqs_df[["protein_id"]],
                    product = seqs_df[["product"]],
                    note = seqs_df[["note"]],
                    seq = as.vector(seqs_df[["translation"]]))
    
    if(sum(lengths(df_list) == 0) > 6) {
      NULL
    } else {
      df_list[lengths(df_list) == 0] <- NA
      
      as.data.frame(df_list)
    }
    
    }, silent = TRUE)
)

save(res, file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiom-results/Gastrointestinal_tract.RData")
# load("Gastrointestinal_tract.RData")
# 
# res_df <- bind_rows(res[sapply(res, function(i) class(i) %in% "data.frame")])
# 
# 
# 
# seqinr::write.fasta(sequences = strsplit(res_df[["seq"]], ""), names = res_df[["protein_id"]], 
#                     file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/raw_seq.fasta")
# 
# 
# files_path <- "/home/michal/Dropbox/signal_final_dat/"
# 
# all_prots <- "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/raw_seq.fasta"
# 
# # SignalP 4.1 --------------------------------------
# 
# signalP_command <- paste0(files_path, "signalp-4.1/signalp -t gram- -f short ", 
#                           all_prots, " > signalP.txt")
# system(signalP_command)
# res_signalP <- read.table("signalP.txt", skip = 1, stringsAsFactors = FALSE) %>% 
#   select(name = V1, 
#          signalP_prob = V9, 
#          signalP_cs = V5) %>% 
#   mutate(signalP_cs = signalP_cs - 1)
# 
# 
# # SignalP 5.0 ---------------------------------------------
# signalP5_command <- paste0(files_path, "signalP5/bin/signalp-5 -organism Gram-negative -params_dir '", 
#                            files_path, "signalP5/bin/model_params' -fasta ", all_prots)
# 
# system(signalP5_command)
# res_signalP5 <- read.table("output_prediction_summary.txt", stringsAsFactors = FALSE, sep = "\t") %>% 
#   select(name = V1, signalP5_prob = V3, signalP5_cs = V7) %>% 
#   mutate(signalP5_prob = as.numeric(signalP5_prob),
#          signalP5_cs = as.character(signalP5_cs),
#          signalP5_cs = sapply(strsplit(signalP5_cs, ".", fixed = TRUE), first),
#          signalP5_cs = sapply(strsplit(signalP5_cs, " "), last),
#          signalP5_cs = sapply(strsplit(signalP5_cs, "-"), first),
#          signalP5_cs = as.numeric(signalP5_cs))
# 
# # DeepSig ----------------------------
# 
# 
# deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ", 
#                           all_prots, " gramn deepSig.txt")
# 
# system(deepSig_command)
# res_deepSig <- read.table("deepSig.txt", stringsAsFactors = FALSE) %>% 
#   mutate(deepSig_prob = ifelse(V2 == "SignalPeptide", V3, 1 - V3),
#          deepSig_cs = as.numeric(ifelse(V4 == "-", NA, V4))) %>% 
#   select(name = V1, 
#          deepSig_prob,
#          deepSig_cs)
# 
# file.remove(c("deepSig.txt", "signalP.txt", "output_prediction_summary.txt"))
# 
# res_sp <- inner_join(res_signalP, res_deepSig) %>% 
#   inner_join(res_signalP5) %>% 
#   rename(protein_id = name) 
# 
# inner_join(res_df, res_sp, by = c("protein_id" = "protein_id")) %>% 
#   write.csv(file = "/home/michal/Dropbox/ann-arbor-collab/gut-microbiome-results/Gastrointestinal_tract.csv", row.names = FALSE)
