args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Argumentos: exp, miembro, init_date (yyyymmddhh)", call.=FALSE)
}

args[1]
args[2]
args[3]

library(metR)
library(data.table)
library(tidyverse)
library(lubridate)
library(here)
source(here("help_functions.R"))
source(here("postprocesamiento.R"))
map_proj <- "+proj=lcc +lat_1=-30.9659996032715 +lat_2=-30.9659996032715 +lat_0=-30.9660034179688 +lon_0=-63.5670013427734 +a=6370000 +b=6370000"

# leo sondeos
files <- list.files(path = "/home/paola.corrales/datosmunin3/DATA/OBS_RELAMPAGO/sondeos_raw",
                    pattern = "cls", full.names = TRUE)

sondeos <- purrr::map(files, ~ read_radiosonde_relampago(.x)) %>%
  rbindlist()

message("Listo sondeos")

# itero sobre pronósticos

files <- list.files(path = paste0("/home/paola.corrales/datosmunin3/EXP/", args[1], "/FCST/", args[3],"/", args[2]), 
                    recursive = TRUE, pattern = "wrfout", full.names = TRUE)

#files <- files[60:70]
for (f in files) {
  
  message(paste("Procesando :", f))
  
  descriptores <- unglue::unglue(f, c("/home/paola.corrales/datosmunin3/EXP/{exp}/FCST/{fcst}/{member}/wrfout_d01_{date}.mean", "/home/paola.corrales/datosmunin3/EXP/{exp}/FCST/{fcst}/{member}/wrfout_d01_{date}"))
#  if (as.numeric(descriptores[[1]][["member"]]) <= 16) {
#	next
#}
  # Leo pronóstico con algo de post procesamiento
  fcst <- ReadNetCDF(f, vars = c(p = "P", "PB", t = "T", qv = "QVAPOR", 
                                 lon = "XLONG", lat = "XLAT")) %>%
    .[, p := p + PB] %>%
    .[, t := tk(t, p)] %>%
    .[, rh := rh(qv, p, t)] %>% 
    .[, td := td(qv, p) + 273.15] %>% 
    .[, ":="(Time = NULL,
             west_east = NULL,
             south_north = NULL,
             qv = NULL,
             PB = NULL)] %>% 
    .[, c("u", "v") := uvmet(f)] %>% 
    .[, ":="(time = ymd_hms(descriptores[[1]][["date"]]),
             init_time = ymd_h(descriptores[[1]][["fcst"]]),
             exp = descriptores[[1]][["exp"]],
             member = descriptores[[1]][["member"]])] %>% 
    .[]
  
  fcst_time <- ymd_hms(descriptores[[1]][["date"]])
  
  intervalo <- interval(fcst_time - minutes(5), fcst_time + minutes(5))
  
  subset <- sondeos[time %within% intervalo] %>% 
    .[, c("xp", "yp") := mesoda::wrf_project(lon, lat, round = FALSE, map_proj = map_proj)]
  
  message(paste0(nrow(subset), " observaciones de sondeos en este tiempo"))
  
  if (nrow(subset) < 1) {
    next
  } 
  
  unique_subset <- subset %>% unique(by = c("lat", "lon"))
  
  unique_site <- unique_subset %>% unique(by = "site")
  
  # intero para cada sondeo
  
  fcst_obs <- purrr::map(unique_site$site, function(x) {
    rx <- range(unique_subset[site == x]$lon, na.rm = TRUE) + c(-1, 1)
    ry <- range(unique_subset[site == x]$lat, na.rm = TRUE) + c(-1, 1)
    out <- fcst %>% 
      .[lat %between% ry & lon %between% rx] %>%
      melt(id.vars = c("bottom_top", "lon", "lat", "time", "init_time", "exp", "member")) %>% 
      .[, c("xp", "yp") := mesoda::wrf_project(lon, lat, map_proj = map_proj)] %>%
      .[, interp_lite(xp, yp, value, 
                      xo = unique_subset[site == x]$xp, 
                      yo = unique_subset[site == x]$yp,
                      output = "points"),
        by = .(bottom_top, variable, time, init_time, exp, member)]
  }) %>% 
    rbindlist() 
  # %>%
  #   .[, c("lon", "lat") := wrf_project(x, y, inverse = TRUE, round = FALSE)]
  message(paste("Interpolación ok!"))
  
  fcst_obs <- fcst_obs[variable != "p"] %>% 
    .[fcst_obs[variable == "p"], on = c("bottom_top", "time", "x", "y", "init_time", "exp", "member")] %>% 
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
    .[, ":="(exp = descriptores[[1]][["exp"]],
             member = descriptores[[1]][["member"]])] %>%
    .[]
  
  fwrite(subset, paste0("/home/paola.corrales/datosmunin3/EXP/", args[1], "/FCST/", args[3], "/sondeos/sondeo_", descriptores[[1]][["exp"]], "_",
                        descriptores[[1]]["member"], "_",  format(fcst_time, "%Y%m%d%H%M%S"), ".csv"))
} 
