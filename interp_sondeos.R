args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Argumentos: exp, miembro, init_date (yyyymmddhh)", call.=FALSE)
}


library(metR)
library(data.table)
library(tidyverse)
library(lubridate)
source("help_functions.R")
source("postprocesamiento.R")

# leo sondeos
files <- list.files(path = "/glade/work/jruiz/sondeos_raw",
                    pattern = "cls", full.names = TRUE)

sondeos <- purrr::map(files, ~ read_radiosonde_relampago(.x)) %>% 
  rbindlist()

message("Listo sondeos")

# itero sobre pronósticos

files <- list.files(path = paste0("/glade/scratch/jruiz/EXP/", args[1], "/FCST/", args[3],"/", args[2]), 
		recursive = TRUE, pattern = "wrfout", full.names = TRUE)
for (f in files) {

message(paste("Procesando :", f))

descriptores <- unglue::unglue(f, c("/glade/scratch/jruiz/EXP/{exp}/FCST/{fcst}/{member}/wrfout_d01_{date}.mean", "/glade/scratch/jruiz/EXP/{exp}/FCST/{fcst}/{member}/wrfout_d01_{date}"))

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

subset <- sondeos[time %within% intervalo]

if (nrow(subset) < 1) {
 	next
} 

unique_subset <- subset %>% unique(by = c("lat", "lon"))

rx <- range(unique_subset$lon, na.rm = TRUE) + c(-1, 1)
ry <- range(unique_subset$lat, na.rm = TRUE) + c(-1, 1)

fcst_obs <- fcst %>% 
  melt(id.vars = c("bottom_top", "lon", "lat", "time", "exp", "member")) %>% 
  .[lat %between% ry & lon %between% rx] %>% 
  .[, interp_lite(lon, lat, value, 
                     xo = unique_subset$lon, yo = unique_subset$lat,
                     output = "points"),
    by = .(bottom_top, variable, time, exp, member)]

#fcst_obs <- fcst %>% melt(id.vars = c("bottom_top", "lon", "lat", "time", "init_time", "exp", "member")) %>% 
#  .[, interp_lite(lon, lat, value, 
#                     xo = unique_subset$lon, yo = unique_subset$lat,
#                     output = "points"), by = .(bottom_top, variable, time, init_time, exp, member)]

fcst_obs <- fcst_obs[variable != "p"] %>% 
  .[fcst_obs[variable == "p"], on = c("bottom_top", "time", "x", "y", "init_time", "exp", "member")] %>% 
  setnames(c("x", "y", "z", "i.z"), c("lon", "lat", "value", "p"))

	
# approx_safe <- purrr::possibly(approx, otherwise = NA)

approx_safe <- function(lon_by, lat_by, variable_by, p) {
sub <- fcst_obs[lon == lon_by & lat == lat_by &
                                   variable == variable_by]
if (nrow(sub) < 2) {
	return(NA_real_)
} else {
	approx(x = sub$p*0.01,  y = sub$value, xout  = p)$y
}

}

subset <- subset %>% 
  melt(measure.vars = c("t", "td", "rh", "u", "v")) %>% 
  .[, fcst_value := approx_safe(.BY$lon, .BY$lat, .BY$variable, p), 
    by = .(lon, lat, variable)] %>%
   .[, ":="(exp = descriptores[[1]][["exp"]],
           member = descriptores[[1]][["member"]])] %>%
   .[]

fwrite(subset, paste0("/glade/scratch/jruiz/EXP/analisis/sondeos/", args[3], "/sondeo_", descriptores[[1]][["exp"]], "_",
			descriptores[[1]]["member"], "_",  format(fcst_time, "%Y%m%d%H%M%S"), ".csv"))
} 
