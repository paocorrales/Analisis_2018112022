# Calculate pp in a latitud band

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)

exp <- "E4"
fcst <- "20181121210000"
files <- Sys.glob(paste0("/glade/scratch/jruiz/EXP/", exp, "/FCST/", fcst, "/wrfout*"))

pp <- lapply(files, function(f) { 
print(f)
a <- ReadNetCDF(f, vars = c("RAINNC", "RAINSH", "RAINC"), subset = list(south_north = 82:130)) %>% 
  .[, .(pp = mean(RAINNC + RAINSH + RAINC, na.rn = TRUE)), by = .(west_east)] %>% 
  .[, exp := exp] %>%
  .[, fcst := fcst] %>%
  .[, date := ymd_hms(str_extract(f, "\\d{4}-\\d{2}-\\d{2}_\\d{2}:\\d{2}:\\d{2}"))]
}) %>%
rbindlist()

fwrite(pp, paste0("pp_hov_fcst", fcst, "_", exp, ".csv"))

