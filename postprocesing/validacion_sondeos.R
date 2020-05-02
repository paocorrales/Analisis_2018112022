library(tidyverse)
library(data.table)
library(metR)

# files <- list.files(path = "/glade/scratch/jruiz/EXP/analisis/sondeos/2018112200/", full.names = TRUE)
files <- list.files(path = "analisis/sondeos/det20181122000000/", full.names = TRUE)

regular_p <- seq(10, 1015, 5)

sondeos <- purrr::map(files, function(f) {

	str <- unglue::unglue(basename(f), "sondeo_{exp}_det{time}.csv")

	out <- fread(f) %>%
		# .[, ens := str[[1]][["member"]]] %>%
		.[, ":="(ele = NULL,
		         azi = NULL,
		         qp = NULL,
		         qt = NULL,
		         qrh = NULL,
		         qu = NULL,
		         qv = NULL,
		         qdZ = NULL)] %>% 
	  setnames(c("value"), c("obs_value")) %>% 
	  .[]

}) %>% rbindlist(fill = TRUE) %>% 
  .[,  fcst_value := if_else(variable %in% c("t", "td"), 
                             fcst_value -273.15, fcst_value)] 
  .[, .(regular_obs = approx(p, obs_value, regular_p)$y,
        regular_fcst = approx(p, fcst_value, regular_p)$y,
        regular_p = regular_p,
        lon = lon[1],
        lat = lat[1]), by = .(launch_time, site, exp, variable)] %>%
  .[,  regular_fcst_mean := mean(regular_fcst, na.rm = TRUE), 
    by = .(variable, regular_p, launch_time, exp)] %>% 
  .[, launch_time := lubridate::ymd_hms(launch_time)]
