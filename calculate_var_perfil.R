# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
library(aiRthermo)


path <- "/glade/scratch/jruiz/EXP/E4/GUESS/2018112*/wrfarw.ensmean"
lat_band <- list(south_north = 105:109)

files <- Sys.glob(path)
levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
levs <- unique(levs)

interpol <- function(x, y, z, ...) {
  IL <- interp::interp(x, y, z, ...)
  CJ(y = IL$y, x = IL$x)[, z := c(IL$z)] 
}

for (f in 1:length(files)) {
  
  date <- lubridate::ymd_hms(str_extract(files[f], "\\d{14}"))
  print(paste("Analisis", date))
  
  ana <- ReadNetCDF(files[f], vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 290)*(P/100000)^(2/7)] %>%
    .[, date := date] 

  ana_sigma <- ana %>% 
    .[, .(T = mean(T, na.rm = TRUE),
          P = mean(P, na.rm = TRUE),
	  QVAPOR = mean(QVAPOR, na.rm = TRUE)), by = .(bottom_top, date)]
  
  ana_lev <- ana %>%
     .[, .(lev = levs, 
           T = approx(P, T, xout = levs*100)$y, 
	   QVAPOR = approx(P, QVAPOR, xout = levs*100)$y), 
         by = .(date, south_north, west_east)] %>%
     .[, .(T = mean(T, na.rm = TRUE),
           QVAPOR = mean(QVAPOR, na.rm = TRUE)), by = .(lev, date)]
 
  if (f == 1) {
    out_sigma <- ana_sigma
    out_lev <- ana_lev
  } else {
    out_sigma <- rbind(out_sigma, ana_sigma)
    out_lev <- rbind(out_lev, ana_lev)
  }

}

fwrite(out_sigma, "/glade/scratch/jruiz/EXP/analisis/perfil_T_Q_guess_E4.csv")
fwrite(out_lev, "/glade/scratch/jruiz/EXP/analisis/perfil_lev_T_Q_guess_E4.csv")

