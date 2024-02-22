library(tidyverse)
library(metR)
library(reticulate)
library(lubridate)

wrf <- import("wrf")
ncdf <- import("netCDF4")
np <- import("numpy")
xr <- import("xarray")



ini_date <- "20181122060000"
exp <- "E6"
fcst_long <- 30
files <- list.files(path = paste0("/glade/scratch/jruiz/EXP/", exp, "/FCST/det", ini_date, "/"), full.names = TRUE)
files <- files[str_detect(files, ":00:00")]

for (f in seq_along(files)[-length(files)]) {

ncfile_ini <- ncdf$Dataset(files[f])
ncfile_end <- ncdf$Dataset(files[f+1])

p <- wrf$getvar(ncfile_ini, "RAINNC")
pp_ini <- wrf$getvar(ncfile_ini, "RAINNC", meta = FALSE) + wrf$getvar(ncfile_ini, "RAINC", meta = FALSE) +
          wrf$getvar(ncfile_ini, "RAINSH", meta = FALSE)
pp_end <- wrf$getvar(ncfile_end, "RAINNC", meta = FALSE) + wrf$getvar(ncfile_end, "RAINC", meta = FALSE) +
          wrf$getvar(ncfile_end, "RAINSH", meta = FALSE)

p$data <- pp_end - pp_ini

date <- as_datetime(ini_date) + hours(f)
dir_out <- paste0("/glade/scratch/jruiz/EXP/analisis/ppacum/", format(as_datetime(ini_date), "%Y%m%d%H"))
path_out <- paste0(dir_out, "/pp_acum_1h_fcst_", exp, "_", as.character(format(date, "%Y%m%d%H%M%S")), ".nc")
#path_out <- paste0("pp_1h_", exp, "_", ini_date, "_f", formatC(f, digits = 2, width = 3, flag = 0), ".nc")
#write_xarray_to_netcdf(p, path_out, engine ="netcdf4")

xarray_array_out <- p$copy(deep = TRUE)
# coordinates are extracted from variable
xarray_array_out$attrs['coordinates'] <- NULL
# wrf-python projection object cannot be processed
xarray_array_out$attrs['projection'] <- as.character(xarray_array_out$attrs['projection'])

xarray_array_out$to_netcdf(path=path_out, mode='w', engine='netcdf4')

ncfile_ini$close()
ncfile_end$close()

message(paste0("Listo: ", files[f+1]))
}



