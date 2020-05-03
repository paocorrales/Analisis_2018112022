library(tidyverse)
library(metR)
library(reticulate)

wrf <- import("wrf")
ncdf <- import("netCDF4")
np <- import("numpy")
xr <- import("xarray")

write_xarray_to_netcdf <- function(xarray_array, output_path, mode = 'w', format = 'NETCDF4', engine = "netcdf4"){
    #writes and xarray in a netcdf format outputfile
    #Uses the xarray typical for wrf-python. The projection objects are transformed into strings
    #to be able to use them as netcdf attributes
    #:param xarray_array: xarray.DataArray
    #:param output_path: str
    #:param format: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT' or 'NETCDF3_CLASSIC'
    #                default: 'NETCDF4'
    #:param group: str, default None
    #:param engine: 'netcdf4', 'scipy' or 'h5netcdf'
    #:param encoding: dict, default: None"

    xarray_array_out <- xarray_array$copy(deep = TRUE)
    # coordinates are extracted from variable
    xarray_array_out$attrs['coordinates'] <- NULL
    # wrf-python projection object cannot be processed
    xarray_array_out$attrs['projection'] <- as.character(xarray_array_out$attrs['projection'])

    xarray_array_out$to_netcdf(path=output_path, mode=mode, format=format, engine=engine)
}

ini_date <- "20181122000000"
exp <- "E4"
files <- list.files(path = paste0("/glade/scratch/jruiz/EXP/", exp, "/FCST/det", ini_date, "/"), full.names = TRUE)

for (f in seq_along(files)-1) {

ncfile_ini <- ncdf$Dataset(files[f])
ncfile_end <- ncdf$Dataset(files[f+1])

p <- wrf$getvar(ncfile_ini, "RAINNC")
pp_ini <- wrf$getvar(ncfile_ini, "RAINNC", meta = FALSE) + wrf$getvar(ncfile_ini, "RAINC", meta = FALSE) + 
          wrf$getvar(ncfile_ini, "RAINSH", meta = FALSE)
pp_end <- wrf$getvar(ncfile_end, "RAINNC", meta = FALSE) + wrf$getvar(ncfile_end, "RAINC", meta = FALSE) + 
          wrf$getvar(ncfile_end, "RAINSH", meta = FALSE) 

p$data <- pp_end - pp_ini

path_out <- paste0("pp_1h_", exp, "_", ini_date, "_f", formatC(f, digits = 2, width = 3, flag = 0), ".nc")
write_xarray_to_netcdf(p, path_out, engine ="netcdf4")

}
