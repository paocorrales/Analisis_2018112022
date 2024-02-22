# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args) != 3) {
#   stop("Argumentos: exp, miembro, init_date (yyyymmddhhmmss)", call.=FALSE)
# }

args <- c("E6_long", "ensmean", "20181108000000")

args[1]
args[2]
args[3]

library(metR)
library(data.table)
library(tidyverse)
library(lubridate)
library(mesoda)
source("/home/paola.corrales/Analisis_2018112022/help_functions.R")
source("/home/paola.corrales/Analisis_2018112022/postprocesamiento.R")

# leo sondeos
files <- list.files(path = "/home/paola.corrales/datosmunin3/DATA/RELAMPAGO/sondeos_raw/",
                    pattern = "cls", full.names = TRUE)

sondeos <- purrr::map(files, ~ read_radiosonde_relampago(.x)) %>%
  rbindlist()

message("Listo sondeos")

# itero sobre pronósticos

files <- Sys.glob(paste0("/home/paola.corrales/datosmunin3/EXP/", args[1], "/ANA/*/analysis.ensmean"))

#files <- files[60:70]
for (f in files) {
  
  message(paste("Procesando :", f))
  
  descriptores <- unglue::unglue(f, c("/home/paola.corrales/datosmunin3/EXP/{exp}/ANA/{date}/analysis.ensmean"))

  fcst_time <- ymd_hms(descriptores[[1]][["date"]])
  
  intervalo <- interval(fcst_time - minutes(30), fcst_time + minutes(30))
  
  subset <- sondeos[time %within% intervalo] %>% 
    .[, c("xp", "yp") := mesoda::wrf_project(lon, lat, round = FALSE)]
  
  message(paste0(nrow(subset), " observaciones de sondeos en este tiempo"))
  
  if (nrow(subset) < 1) {
    next
  } 
  
  # Leo pronóstico con algo de post procesamiento
  fcst <- ReadNetCDF(f, vars = c(p = "P", "PB", t = "T", qv = "QVAPOR",
                                 lon = "XLONG", lat = "XLAT")) %>%
    .[, p := p + PB] %>%
    .[, t := tk(t, p, T_BASE = 300)] %>%
    .[, rh := rh(qv, p, t)] %>%
    .[, td := td(qv, p) + 273.15] %>%
    .[, ":="(Time = NULL,
             west_east = NULL,
             south_north = NULL,
             qv = NULL,
             PB = NULL)] %>%
    .[, c("u", "v") := uvmet(f)] %>%
    .[, ":="(date = ymd_hms(descriptores[[1]][["date"]]),
             exp = descriptores[[1]][["exp"]])] %>%
    .[]


  unique_subset <- subset %>% unique(by = c("lat", "lon"))
  
  unique_site <- unique_subset %>% unique(by = "site")
  
  # intero para cada sondeo
  
  fcst_obs <- purrr::map(unique_site$site, function(x) {
    rx <- range(unique_subset[site == x]$lon, na.rm = TRUE) + c(-1, 1)
    ry <- range(unique_subset[site == x]$lat, na.rm = TRUE) + c(-1, 1)
    out <- fcst %>% 
      .[lat %between% ry & lon %between% rx] %>%
      melt(id.vars = c("bottom_top", "lon", "lat", "date", "exp")) %>% 
      .[, c("xp", "yp") := mesoda::wrf_project(lon, lat)] %>%
      .[, interp_lite(xp, yp, value, 
                      xo = unique_subset[site == x]$xp, 
                      yo = unique_subset[site == x]$yp,
                      output = "points"),
        by = .(bottom_top, variable, date, exp)]
  }) %>% 
    rbindlist() 
  # %>%
  #   .[, c("lon", "lat") := wrf_project(x, y, inverse = TRUE, round = FALSE)]
  message(paste("Interpolación ok!"))
  
  fcst_obs <- fcst_obs[variable != "p"] %>% 
    .[fcst_obs[variable == "p"], on = c("bottom_top", "date", "x", "y", "exp")] %>% 
    setnames(c("x", "y", "z", "i.z"), c("xp", "yp", "value", "p"))
  
  approx_safe <- function(lon_by, lat_by, variable_by, p) {
    sub <- fcst_obs[xp == lon_by & yp == lat_by &
                      variable == variable_by]
    if (nrow(sub) < 2) {
      return(NA_real_)
    } else {
      approx(x = sub$p*0.01,  y = sub$value, xout  = p)$y
    }
    
  }
  
  subset <- subset %>% 
    melt(measure.vars = c("t", "td", "rh", "u", "v")) %>% 
    .[, fcst_value := approx_safe(.BY$xp, .BY$yp, .BY$variable, p), 
      by = .(xp, yp, variable)] %>%
    .[, ":="(exp = descriptores[[1]][["exp"]])] %>%
    .[]
  
  fwrite(subset, paste0("/home/paola.corrales/datosmunin3/EXP/derived_data/sondeos/ANA/sondeo_", descriptores[[1]][["exp"]], "_",
              format(fcst_time, "%Y%m%d%H%M%S"), "new.csv"))
} 
