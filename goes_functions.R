
# projection --------------------------------------------------------------

goes_projection <- function(x, y, ncfile) {

    proj_info <- ncatt_get(ncfile, "goes_imager_projection")
  h <- proj_info$perspective_point_height 
  
  map_proj <- paste0("+proj=geos",
                     " +h=",  proj_info$perspective_point_height,
                     " +lon_0=", proj_info$longitude_of_projection_origin,
                     " +sweep=", proj_info$sweep_angle_axis,
                     " +ellps=GRS80")
  
  x_atr <- ncatt_get(ncfile, "x")
  x <- (x*x_atr$scale_factor + x_atr$add_offset)*h
  y_atr <- ncatt_get(ncfile, "y")
  y <- (y*y_atr$scale_factor + y_atr$add_offset)*h

  proj4::project(list(x, y), map_proj, inverse = TRUE)
}

# ggplot2::ggplot() +
#   
#   scattermore::geom_scattermost(data_sub[!is.na(z_norm), .(lon, lat)], interpolate = TRUE, pixels = c(512, 512),
#                                 viridisLite::viridis(100)[1 + 99*data_sub[!is.na(z_norm), z_norm]]) +
#   geom_path(data = map, aes(long, lat, group = group)) +
#   ggplot2::coord_equal(ylim = c(-45, -20), xlim = c(-77, -40))


# Planck ------------------------------------------------------------------

rad_to_tb <- function(rad, mu, h = 6.629e-34, c = 2.998e8, kb = 1.381e-23) {
  mu <- c(mu)
  rad <- c(rad)
  tb <- h*c/(kb*mu*log(1 + 2*h*c^2/(rad*mu^5)))
  
  return(tb)
}

