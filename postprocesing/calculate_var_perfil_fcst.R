# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)
# source("~/Analisis_2018112022/postprocesamiento.R")

wrf_path <- "/home/paola.corrales/datosmunin3/EXP/"

lon_band = c(-66.5, -61.5)
lat_band = c(-35.5, -29)

exp <- "E6"
run <- "fcst"
ini_date <- ymd_hms("20181122000000")
ciclos <- 36

# dates <- seq.POSIXt(ini_date, by = "hour",
#                     length.out = ciclos)

levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
levs <- unique(levs)

rh <- function(QVAPOR, T, P) {
  
  # T en kelvin, P en Pa
  # QVAPOR dimensionless
  
  rh <- (0.263 * P * QVAPOR) / (exp((17.67*(T - 273.16))/(T - 29.65)))
  return(rh)
}

lead_time <- seq(0, ciclos)

for (exp in exps) {
  wrf <- purrr::map(lead_time, function(l) {
    
    print(paste("Pronostico", l))
    
    date <- ini_date + hours(l)
    file_fcst <- Sys.glob(paste0(fcst_path, "/", exp, "/FCST/", format(ini_date, "%Y%m%d%H"),
                                 "_det/wrfout_d02_",
                                 format(date, "%Y"), "-",
                                 format(date, "%m"), "-",
                                 format(date, "%d"), "_",
                                 format(date, "%H"), ":",
                                 format(date, "%M"), ":",
                                 format(date, "%S")))

    
    if (!file.exists(file_fcst)) {
      return(NULL)
    }
    
    ReadNetCDF(file_fcst, vars = c(p = "P", "PB", t = "T", q = "QVAPOR", 
                                   lon = "XLONG", lat = "XLAT")) %>%
      .[, p := p + PB] %>%
      .[, t := tk(t, p)] %>%
      .[, ":="(Time = NULL,
               PB = NULL, 
               date = date)] %>% 
      .[, c("u", "v") := uvmet(file_fcst)] %>% 
      .[lat %between% lat_band & lon %between% lon_band] %>% 
      .[, .(lev = levs, 
            t = approx(p, t, xout = levs*100)$y, 
            q = approx(p, q, xout = levs*100)$y,
            u = approx(p, u, xout = levs*100)$y,
            v = approx(p, v, xout = levs*100)$y), 
        by = .(date, south_north, west_east)] %>% 
      .[, .(t = mean(t, na.rm = TRUE),
            q = mean(q, na.rm = TRUE),
            u = mean(u, na.rm = TRUE),
            v = mean(v, na.rm = TRUE)), by = .(date, lev)] %>% 
      .[, exp := exp]
    
    
  }) %>%
    rbindlist()
  
  
  fwrite(wrf, paste0(wrf_path, "/derived_data/perfiles_fcst_det_", exp, "_", format(ini_date, "%Y%m%d%H%M%S"), ".csv"))
  
}
