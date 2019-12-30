library(tidyverse)
library(data.table)

string <- read_lines("prepobs_errtable.global")

split_chunks <- function(string, pattern = "OBSERVATION TYPE", negate = FALSE) {
  start <- stringr::str_which(string, pattern, negate = negate)
  end <- c(start-1, length(string))[-1]
  
  obs_type <- stringr::str_extract(string[start], "\\d+")
  
  chunks <- purrr::map(seq_along(start), ~ string[(start[.x]+1):end[.x]])
  
  names(chunks) <- obs_type  
  chunks
}

string <- split_chunks(errtable, "OBSERVATION TYPE")

chunk <- string[[1]]

parse_chunk <- function(chunk) {
  str_split(chunk, " ") %>% 
    map( ~ str_subset(.x, "^$", negate = TRUE) %>% 
           setNames(c("P", "T", "Q", "U", "V", "PW")) %>%
           as.numeric()) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
}

errtable <- map(string, parse_chunk) %>% 
  rbindlist(idcol = "type") 

str(errtable)
