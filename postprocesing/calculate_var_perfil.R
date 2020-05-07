# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
source("../postprocesamiento.R")

path <- "/glade/scratch/jruiz/EXP/E6/ANA/2018112*/analysis.ensmean"
lat_band <- list(south_north = 40:299)
lon_band <- list(west_east = 70:180)

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
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, c("U", "V") := uvmet(files[f])] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>% 
    .[, date := date] 

  #ana_sigma <- ana %>% 
  #  .[, .(T = mean(T, na.rm = TRUE),
  #        P = mean(P, na.rm = TRUE),
  #	  QVAPOR = mean(QVAPOR, na.rm = TRUE)), by = .(bottom_top, date)]
  
  ana_lev <- ana %>%
     .[, .(lev = levs, 
           T = approx(P, T, xout = levs*100)$y, 
	   QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
           U = approx(P, U, xout = levs*100)$y,
           V = approx(P, V, xout = levs*100)$y,
           SPD = approx(P, SPD, xout = levs*100)$y), 
         by = .(date, south_north, west_east)] %>%
     .[, .(T = mean(T, na.rm = TRUE),
           QVAPOR = mean(QVAPOR, na.rm = TRUE),
           U = mean(U, na.rm = TRUE),
           V = mean(V, na.rm = TRUE),
           SPD = mean(SPD, na.rm = TRUE)), by = .(lev, date)]
 
  if (f == 1) {
  #  out_sigma <- ana_sigma
    out_lev <- ana_lev
  } else {
  #  out_sigma <- rbind(out_sigma, ana_sigma)
    out_lev <- rbind(out_lev, ana_lev)
  }

}

#fwrite(out_sigma, "/glade/scratch/jruiz/EXP/analisis/perfil_T_Q_ana_E4.csv")
fwrite(out_lev, "/glade/scratch/jruiz/EXP/analisis/perfiles_ana_E6.csv")

