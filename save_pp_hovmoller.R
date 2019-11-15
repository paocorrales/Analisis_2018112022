# Calculate pp in a latitud band

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)

a <- ReadNetCDF("E3/ANA/20181122060000/analysis.ensmean", vars = c("RAINNC", "RAINSH", "RAINC"), subset = list(south_north = 82:130)) %>% 
  .[, .(pp = mean(RAINNC + RAINSH + RAINC, na.rn = TRUE)), by = .(west_east)] %>% 
  .[]
a

