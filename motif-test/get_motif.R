library(dplyr)

aas <- toupper(biogram:::return_elements("prot"))




process_single_motif <- function(motif, len_min, len_max) {
  raw_motif1 <- gsub("[[:space:]]", "", motif) %>% 
    gsub("X", paste0(toupper(biogram:::return_elements("prot")), collapse = ""), .) %>% 
    strsplit("") %>% 
    unlist
  
  motif <- if("-" %in% raw_motif1) {
    
    lhs <- raw_motif1[1L:(which(raw_motif1 == "-") - 1)]
    rhs <- raw_motif1[(which(raw_motif1 == "-") +1):length(raw_motif1)]
    
    setdiff(lhs, rhs)
  } else {
    raw_motif1
  }
  
  paste0("[", paste0(motif, collapse = ""), "]{", len_min, ",", len_max, "}")
}



grep(pattern = paste0("B", process_single_motif("A", 2, 3), "B"), x = c("BAB", "BAAB", "BAAAB", "BAAAAAB"))

     