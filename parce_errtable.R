library(tidyverse)
library(data.table)

errtable <- read_lines("prepobs_errtable.global")

errtable %>% 
  str_which("OBSERVATION TYPE")

errtable %>%  
  str_split(" ") %>% 
  map(~ str_subset(.x, "^$", negate = TRUE)) %>% 
  
