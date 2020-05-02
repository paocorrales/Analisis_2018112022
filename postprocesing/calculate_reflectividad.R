library(reticulate)
wrf <- import("wrf")
ncdf <- import("netCDF4")

get_dbz <- function(file_in, file_out) {
  
  ncfile <-  ncdf$Dataset(file_in)
  dbz <- wrf$g_dbz$get_max_dbz(ncfile)
  
  xarray_array_out <-  dbz$copy(deep = "True")
  xarray_array_out$attrs['coordinates'] <- NULL
  
  xarray_array_out$attrs['projection'] = as.character(xarray_array_out$attrs['projection']$projection)
  
  xarray_array_out$to_netcdf(path = file_out, 
                             mode = "w",
                             engine = "netcdf4")
  return(file_out)
}


files <- list.files("/glade/scratch/jruiz/EXP/E3/ANA", recursive = TRUE, full.names = TRUE, pattern = "analysis.ensmean")

for (file in files) {
  date <- basename(dirname(file))
  exp <- basename(dirname(dirname(dirname(file))))
  out_file <- file.path("dbz", paste0("maxdbz_ana_", exp, "_", date, ".nc"))
  
  message("Proccessing", file)
  
  get_dbz(file, out_file)
}
