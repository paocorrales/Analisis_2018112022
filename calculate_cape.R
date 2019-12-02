# Calculate CAPE in a latitud band

library(tidyverse)
library(data.table)
library(metR)
library(aiRthermo)

capear <- function(p, t, q) {
  # browser()
  cape <- CAPE_CIN(p, t, q)
  
  list(cape = cape$cape,
       cin = cape$cin,
       outCode = cape$outCode)
}

path <- "/glade/scratch/jruiz/EXP/E4/ANA/2018112*/analysis.ensmean"
lat_band <- list(south_north = 105:109)

files <- Sys.glob(path)
files

for (f in 1:length(files)) {
  
  date <- lubridate::ymd_hms(str_extract(files[f], "\\d{14}"))
  print(paste("Analisis", date))
  
  ana <- ReadNetCDF(files[f], vars = c("P", "PB", "T", "QVAPOR"), subset = lat_band) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 290)*(P/100000)^(2/7)] %>%
    .[, date := date] %>% 
    .[, c("PB") := NULL]
  
  invisible(capture.output(cape <- ana[, capear(P, T, QVAPOR), by = .(south_north, west_east, date)]))
  
  if (f == 1) {
    out <- cape
  } else {
    out <- rbind(out, cape)
  }
}

fwrite(out, "/glade/scratch/jruiz/EXP/analisis/cape_cin_ana_E4.csv")
