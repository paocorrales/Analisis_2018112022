# Calculate RMSD vertical profile 

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)
source("../postprocesamiento.R")

path <- "/home/paola.corrales/datosmunin3/EXP/"

lat_band <- c(40, 299)
lon_band <- c(70, 180)

exp <- "E2"
ini_date <- ymd_hms("20181120180000")
ciclos <- 67

dates <- seq.POSIXt(ini_date, by = "hour",
                    length.out = ciclos)

levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
levs <- unique(levs)

out_rmsd <- list()

for (d in 1:length(dates)) {
  
  # date <- lubridate::ymd_hms(str_extract(files[f], "\\d{14}"))
  
  print(paste("Analisis", dates[d]))
  
  ana_file <- paste0(path, exp, "/ANA/", format(dates[d], "%Y%m%d%H%M%S"), "/analysis.ensmean")
  
  ana <- ReadNetCDF(ana_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, c("U", "V") := uvmet(ana_file)] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>% 
    .[, date := dates[d]] %>% 
    .[south_north %between% lat_band & west_east %between% lon_band] %>% 
    .[, .(lev = levs, 
          T = approx(P, T, xout = levs*100)$y, 
          QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
          U = approx(P, U, xout = levs*100)$y,
          V = approx(P, V, xout = levs*100)$y,
          SPD = approx(P, SPD, xout = levs*100)$y), 
      by = .(date, south_north, west_east)] %>% 
    melt(id.vars = c("date", "south_north", "west_east", "lev"),
         value.name = "ana")
  
  
  guess_file <- paste0(path, exp, "/GUESS/", format(dates[d], "%Y%m%d%H%M%S"), "/wrfarw.ensmean")
  
  guess <- ReadNetCDF(guess_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, c("U", "V") := uvmet(guess_file)] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>% 
    .[, date := dates[d]] %>% 
    .[south_north %between% lat_band & west_east %between% lon_band] %>% 
    .[, .(lev = levs, 
          T = approx(P, T, xout = levs*100)$y, 
          QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
          U = approx(P, U, xout = levs*100)$y,
          V = approx(P, V, xout = levs*100)$y,
          SPD = approx(P, SPD, xout = levs*100)$y), 
      by = .(date, south_north, west_east)] %>% 
    melt(id.vars = c("date", "south_north", "west_east", "lev"),
         value.name = "guess")
  
  
 mse <- cbind(ana, guess) %>% 
  .[, .(rmsd = sqrt(mean((ana - guess)^2, na.rm = TRUE))), by = .(lev, variable, date)]
  
  out_rmsd[[d]] <- mse
  
}


fwrite(rbindlist(out_rmsd), paste0(path, "/derived_data/perfiles_rmsd_", exp, ".csv"))

