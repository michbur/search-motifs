library(dplyr)
library(genbankr)
library(Biostrings)
library(pbapply)

#awk -v n=1 '/^$/{close("out"n);n++;next} {print > "/home/michal/Dropbox/ann-arbor-collab/gut-microbiom-data/out"n}' 42773.ADKP01000001-ADKP01000236.nuc.gbk

res <- pblapply(list.files("/home/michal/Dropbox/ann-arbor-collab/gut-microbiom-data/", full.names = TRUE)[1L:3], function(file_name) 
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

save(res, file = "Gastrointestinal_tract.RData")