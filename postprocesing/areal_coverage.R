library(metR)
library(tidyverse)
library(data.table)
library(lubridate)
library(unglue)
library(mesoda)
# source("/glade/scratch/pcorrales/Analisis_2018112022/help_functions.R")
# source("/glade/scratch/pcorrales/Analisis_2018112022/postprocesamiento.R")


wrf_path <- "~/datosmunin3/EXP/"
out_path <- "~/datosmunin3/EXP/derived_data/ppacum/"
exp <- "EG3_long"
run <- "ana"

ini_date <- ymd_hms("20181108000000")
ciclos <- 379

acumulado <- 1
q <- c(1, 5, 10, 25, 50)

lon_band <- c(-66.5, -61.5)
lat_band <- c(-35.5, -29)


dates <- seq.POSIXt(ini_date + hours(acumulado), by = "hour",
                    length.out = ciclos - acumulado)

pp_wrf_out <- purrr::map_dbl(dates, function(d) {
  
  date <- d
  
  print(d)
  
  
  files_wrf <-  purrr::map(seq(0, acumulado - 1), function(l) {
    list.files(path = paste0(wrf_path, exp, "/", toupper(run), "/", format(date - hours(l), "%Y%m%d%H%M%S")),
               full.names = TRUE,
               recursive = TRUE,
               pattern = "analysis.mem*")
  }) %>% unlist()
  
  
  pp_wrf <- purrr::map(files_wrf, function(f) {
    
    metadatos <- unglue(basename(f), "analysis.mem0{mem}")
    
    ReadNetCDF(f, vars = c("RAINNC", "RAINC", "RAINSH",
                           lon = "XLONG", lat = "XLAT")) %>%
      .[, ":="(pp_acum = RAINNC + RAINC + RAINSH,
               exp = exp,
               date = date,
               mem = metadatos[[1]][["mem"]])] %>%
      .[, ":="(RAINNC = NULL,
               RAINC = NULL,
               RAINSH = NULL,
               Time = NULL)] %>%
      .[lon %between% lon_band & lat %between% lat_band]
  }) %>%
    rbindlist() %>%
    .[, .(pp_acum = sum(pp_acum)), by = .(exp, mem, date, lon, lat)] %>%
    .[, c("x", "y") := wrf_project(lon, lat)]
  
  pp_prop <- purrr::map(q, function(f){
    
    pp_wrf %>%
      .[, .(area = sum(pp_acum > f)/.N), by = .(mem, exp, date)] %>%
      .[, umbral := f]
  }) %>%
    rbindlist()
  
  
  saveRDS(pp_prop, paste0(out_path, exp, "_ana_", format(date, "%Y%m%d%H"), "_area_acum_", acumulado, "h_new.rds"))
  
  d
})

  