# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)
source("~/Analisis_2018112022/postprocesamiento.R")

path <- "/home/paola.corrales/datosmunin3/EXP/"

lat_band <- c(40, 299)
lon_band <- c(70, 180)

exp <- "E9"
ini_date <- ymd_hms("20181120180000")
ciclos <- 67

dates <- seq.POSIXt(ini_date, by = "hour",
                    length.out = ciclos)

levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
levs <- unique(levs)

rh <- function(QVAPOR, T, P) {
  
  # T en kelvin, P en Pa
  # QVAPOR dimensionless
  
  rh <- (0.263 * P * QVAPOR) / (exp((17.67*(T - 273.16))/(T - 29.65)))
 return(rh)
}

for (d in 1:length(dates)) {
  print(paste("Analisis", dates[d]))
  
  ana_file <- paste0(path, exp, "/ANA/", format(dates[d], "%Y%m%d%H%M%S"), "/analysis.ensmean")
  
  ana <- ReadNetCDF(ana_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, RH := rh(QVAPOR, T, P)] %>% 
    .[, c("U", "V") := uvmet(ana_file)] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>% 
    .[, date := dates[d]] %>% 
    .[south_north %between% lat_band & west_east %between% lon_band] %>% 
    .[, .(lev = levs, 
          T = approx(P, T, xout = levs*100)$y, 
          RH = approx(P, RH, xout = levs*100)$y,
          QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
          U = approx(P, U, xout = levs*100)$y,
          V = approx(P, V, xout = levs*100)$y,
          SPD = approx(P, SPD, xout = levs*100)$y), 
      by = .(date, south_north, west_east)] %>%
     .[, .(T = mean(T, na.rm = TRUE),
           RH = mean(RH, na.rm = TRUE),
           QVAPOR = mean(QVAPOR, na.rm = TRUE),
           U = mean(U, na.rm = TRUE),
           V = mean(V, na.rm = TRUE),
           SPD = mean(SPD, na.rm = TRUE)), by = .(lev, date)]
  
  # guess_file <- paste0(path, exp, "/GUESS/", format(dates[d], "%Y%m%d%H%M%S"), "/wrfarw.ensmean")
  # 
  # guess <- ReadNetCDF(guess_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
  #   .[, P := P + PB] %>%
  #   .[, T := (T + 300)*(P/100000)^(2/7)] %>%
  #   .[, c("U", "V") := uvmet(ana_file)] %>%
  #   .[, SPD := sqrt(U^2 + V^2)] %>% 
  #   .[, date := dates[d]] %>% 
  #   .[south_north %between% lat_band & west_east %between% lon_band] %>% 
  #   .[, .(lev = levs, 
  #         T = approx(P, T, xout = levs*100)$y, 
  #         QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
  #         U = approx(P, U, xout = levs*100)$y,
  #         V = approx(P, V, xout = levs*100)$y,
  #         SPD = approx(P, SPD, xout = levs*100)$y), 
  #     by = .(date, south_north, west_east)] %>%
  #   .[, .(T = mean(T, na.rm = TRUE),
  #         QVAPOR = mean(QVAPOR, na.rm = TRUE),
  #         U = mean(U, na.rm = TRUE),
  #         V = mean(V, na.rm = TRUE),
  #         SPD = mean(SPD, na.rm = TRUE)), by = .(lev, date)]
  # 
  if (d == 1) {
    # out_guess <- guess
    out_ana <- ana
  } else {
    # out_guess <- rbind(out_guess, guess)
    out_ana <- rbind(out_ana, ana)
  }

}

# fwrite(out_guess, paste0(path, "/derived_data/perfiles_guess_", exp, ".csv"))
fwrite(out_ana, paste0(path, "/derived_data/perfiles_ana_", exp, "_new.csv"))

