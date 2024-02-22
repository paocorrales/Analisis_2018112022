# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args) != 2) {
#   stop("Argumentos: exp, init_date (yyyymmddhhmmss)", call.=FALSE)
# }

args <- c("E9_long", "20181122000000")



library(metR)
library(data.table)
library(tidyverse)
library(lubridate)
source("~/Analisis_2018112022//help_functions.R")
source("~/Analisis_2018112022//postprocesamiento.R")

# leo sondeos
files <- list.files(path = "/home/paola.corrales/datosmunin3/DATA/RELAMPAGO/sondeos_raw/",
                    pattern = "cls", full.names = TRUE)

sondeos <- purrr::map(files, ~ read_radiosonde_relampago(.x)) %>%
  rbindlist()

message("Listo sondeos")

# itero sobre pronósticos

first_date <- ymd_hms("20181108120000")
fcsts <- 16
# lead_time <- 36

# plan(multisession, workers = 18)

f = 0
while (f <= fcsts) {
  
  ini_date <- first_date + days(f)
  message(paste0("Arranca ", ini_date))

files <- Sys.glob(paste0("~/datosmunin3/EXP/", args[1], "/FCST/", format(ini_date, "%Y%m%d%H%M%S"),"/wrfout_d02*"))

out <- list()
for (fi in 1:length(files)) {

  message(paste("Procesando :", files[fi]))
  
  descriptores <- unglue::unglue(files[fi], c("/home/paola.corrales/datosmunin3/EXP/{exp}/FCST/{fcst}/wrfout_d02_{date}", "/glade/scratch/jruiz/EXP/{exp}/FCST/det{fcst}/wrfout_d01_{date}"))
  
  fcst_time <- ymd_hms(descriptores[[1]][["date"]])
  
  intervalo <- interval(fcst_time - minutes(30), fcst_time + minutes(30))
  
  subset <- sondeos[time %within% intervalo & lon %between% c(-67.68781, -53.76321)] %>% 
    .[, c("xp", "yp") := mesoda::wrf_project(lon, lat, round = FALSE)]
  
  message(paste0(nrow(subset), " observaciones de sondeos en este tiempo"))
  
  if (nrow(subset) < 1) {
    next
  } 
  
  
# Leo pronóstico con algo de post procesamiento
  fcst <- ReadNetCDF(files[fi], vars = c(p = "P", "PB", t = "T", qv = "QVAPOR", 
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
    .[, c("u", "v") := uvmet(files[fi])] %>% 
    .[, ":="(time = ymd_hms(descriptores[[1]][["date"]]),
             init_time = ymd_hms(descriptores[[1]][["fcst"]]),
             exp = descriptores[[1]][["exp"]])] %>%
    .[]
  
  
  unique_subset <- subset %>% unique(by = c("lat", "lon"))
  
  unique_site <- unique_subset %>% unique(by = "site") %>% 
    .[lon %between% c(-67.5, -53.5) & lat %between% range(fcst$lat)]
  
  # intero para cada sondeo
  
  fcst_obs <- purrr::map(unique_site$site, function(x) {
    message(x)
    rx <- range(unique_subset[site == x]$lon, na.rm = TRUE) + c(-0.06, 0.06)
    ry <- range(unique_subset[site == x]$lat, na.rm = TRUE) + c(-0.06, 0.06)
    x_out <-  unique_subset[site == x]$xp
    y_out <-  unique_subset[site == x]$yp
    
    out <- fcst %>% 
      .[lat %between% ry & lon %between% rx] %>%
      melt(id.vars = c("bottom_top", "lon", "lat", "time", "init_time", "exp")) %>% 
      .[, c("xp", "yp") := mesoda::wrf_project(lon, lat, round = c(-1, -1))] %>%
      .[, metR::Interpolate(value ~ xp + yp, 
                            x.out = x_out, 
                            y.out = y_out,
                            grid = FALSE),
        by = .(bottom_top, variable, time, init_time, exp)] 
  }) %>% 
    rbindlist() %>% 
    setnames(c("xp", "yp", "value"), c("x", "y", "z"))
  # %>%
  #   .[, c("lon", "lat") := wrf_project(x, y, inverse = TRUE, round = FALSE)]
  message(paste("Interpolación ok!"))
  
  fcst_obs <- fcst_obs[variable != "p"] %>% 
    .[fcst_obs[variable == "p"], on = c("bottom_top", "time", "x", "y", "init_time", "exp")] %>% 
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
  
  out[[fi]] <- subset %>% 
    melt(measure.vars = c("t", "td", "rh", "u", "v")) %>% 
    .[, fcst_value := approx_safe(.BY$xp, .BY$yp, .BY$variable, p), 
      by = .(xp, yp, variable)] %>%
    .[, ":="(exp = descriptores[[1]][["exp"]])]

}

rbindlist(out) %>% 
  fwrite(paste0("~/datosmunin3/EXP/derived_data/sondeos/sondeo_", descriptores[[1]][["exp"]], "_",
                      "fcst_d02_",  format(ini_date, "%Y%m%d%H%M%S"), ".csv"))

f = f + 1
}
