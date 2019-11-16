library(reticulate)
wrf <- import("wrf")
ncdf <- import("netCDF4")

dale_cape <- function(file_in, file_out) {
  
  ncfile <-  ncdf$Dataset(file_in)
  cape <- wrf$g_cape$get_2dcape(ncfile)
  
  xarray_array_out <-  cape$copy(deep = "True")
  xarray_array_out$attrs['coordinates'] <- NULL
  
  xarray_array_out$attrs['projection'] = as.character(xarray_array_out$attrs['projection']$projection)
  
  xarray_array_out$to_netcdf(path = file_out, 
                             mode = "w",
                             engine = "netcdf4")
  return(file_out)
}


files <- list.files("/glade/scratch/jruiz/EXP/E4/ANA", recursive = TRUE, full.names = TRUE, pattern = "analysis.ensmean")

for (file in files) {
  date <- basename(dirname(file))
  exp <- basename(dirname(dirname(dirname(file))))
  out_file <- file.path("mucape", paste0("mcape_ana_", exp, "_", date, ".nc"))
  
  message("Proccessing", file)
  
  dale_cape(file, out_file)
}
