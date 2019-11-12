# Calculate CAPE in a latitud band

library(tidyverse)
library(data.table)
library(metR)
library(aiRthermo)


path <- "/glade/scratch/jruiz/EXP/E3/ANA/2018112*/analysis.ensmean"
lat_band <- list(south_north = 105:109)

files <- Sys.glob(path)

for (f in 1:length(files)) {
  
  date <- lubridate::ymd_hms(str_extract(files[f], "\\d{14}"))
  print(paste("Analisis", date))
  
  ana <- ReadNetCDF(files[f], vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 290)*(P/100000)^(2/7)] %>%
    .[, date := date] %>% 
    .[, .(T = mean(T, na.rm = TRUE),
		  QVAPOR = mean(QVAPOR, na.rm = TRUE)), by = .(bottom_top, date)]
  
  if (f == 1) {
    out <- ana
  } else {
    out <- rbind(out, ana)
  }
}

fwrite(out, "/glade/scratch/jruiz/EXP/analisis/perfil_T_Q_ana_E3.csv")
