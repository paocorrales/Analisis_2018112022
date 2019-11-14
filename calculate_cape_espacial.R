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

fcst <- c("20181121210000", "20181122000000", "20181122030000", "20181122060000")

files <- map(fcst, function(x) paste0("/glade/scratch/jruiz/EXP/E3/ANA/", x, "/analysis.ensmean"))
lat_band <- list(south_north = 105:109)

files

for (f in 1:length(files)) {
  
  date <- lubridate::ymd_hms(str_extract(files[f], "\\d{14}"))
  print(paste("Analisis", date))
  
  ana <- ReadNetCDF(files[[f]], vars = c("P", "PB", "T", "QVAPOR")) %>%
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

fwrite(out, "/glade/scratch/jruiz/EXP/analisis/cape_cin_espacial_ana_E3.csv")
