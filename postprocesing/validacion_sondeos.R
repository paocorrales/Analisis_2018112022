library(tidyverse)
library(data.table)
library(metR)

files <- list.files(path = "/glade/scratch/jruiz/EXP/analisis/sondeos/2018112200/", full.names = TRUE)

sondeos <- purrr::map(files, function(f) {

	str <- unglue::unglue(basename(f), "sondeo_{exp}_{member}_{time}.csv")

	out <- fread(f) %>%
		.[, ens := str[[1]][["member"]]] %>%
		.[]

}) %>% rbindlist(fill = TRUE)

head(sondeos)

