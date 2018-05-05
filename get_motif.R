library(dplyr)

aas <- LETTERS[1L:5]

raw_motif1 <- gsub("[[:space:]]", "", "X - BC") %>% 
  strsplit("") %>% 
  unlist

#if("-" %in% raw_motif1)

lhs <- raw_motif1[1L:(which(raw_motif1 == "-") - 1)]
rhs <- raw_motif1[(which(raw_motif1 == "-") +1):length(raw_motif1)]

if(grepl("X", lhs)) {
  lhs <- aas
}

setdiff(lhs, rhs)
