library(dplyr)

aas <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
         "P", "Q", "R", "S", "T", "V", "W", "Y")


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
