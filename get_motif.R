library(dplyr)

aas <- LETTERS[1L:5]

raw_motif1 <- "X - B"
gsub("[[:space:]]", "", raw_motif1) %>% 
  strsplit("") %>% 
  unlist


