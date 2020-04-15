
# projection --------------------------------------------------------------

goes_project <- function(lon, lat, inverse = FALSE, round = TRUE) {
  # map_proj <- "+proj=geos +h=35785831.0 +lon_0=-75 +sweep=x +R=6378137"
  map_proj <- "+proj=geos +a=6378137. +b=6356752.31414 +lon_0=-89.5 +f=.00335281068119356027 +h=35786023. +sweep=x"
  xy <- proj4::project(list(lon, lat), map_proj, inverse = inverse)
  if (round == TRUE) {
    list(x = round(xy$x, -3),
         y = round(xy$y, -4))
  } else {
    list(x = xy$x,
         y = xy$y)
  }
}

library(ncdf4)
file <- "~/Downloads/OR_ABI-L1b-RadF-M3C07_G16_s20183261200367_e20183261211145_c20183261211181 (1).nc"
ncfile <- nc_open(file)

z <- ncvar_get(ncfile, "Rad")

proj_info <- ncatt_get(ncfile, "goes_imager_projection")
h <- proj_info$perspective_point_height 

map_proj <- paste0("+proj=geos",
                   " +h=",  proj_info$perspective_point_height,
                   " +lon_0=", proj_info$longitude_of_projection_origin,
                   " +sweep=", proj_info$sweep_angle_axis,
                   " +ellps=GRS80")


# Data info
x_atr <- ncatt_get(ncfile, "x")
x <- ncvar_get(ncfile, "x")*x_atr$scale_factor + x_atr$add_offset
y_atr <- ncatt_get(ncfile, "y")
y <- ncvar_get(ncfile, "y")*y_atr$scale_factor + y_atr$add_offset


library(data.table)
data <- data.table::CJ(y, x, sorted = FALSE)[, c(2, 1)]


data[, z := c(z)]
data[z > 1, z := NA]
data[, c("x_m", "y_m") := list(x*h, y*h) ]




library(data.table)
library(ggplot2)
map <- fortify(rnaturalearth::ne_countries())
map <- as.data.table(map)
map[, c("x", "y") := proj4::project(list(long, lat), map_proj, ellps.default = "GRS80")]

data[, z_norm := scales::rescale(z, to = c(0, 1))]
ggplot2::ggplot() +
  scattermore::geom_scattermost(data[!is.na(z_norm), .(x_m, y_m)],
                                viridisLite::viridis(100)[1 + 99*data[!is.na(z_norm), z_norm]]) +
  # geom_hline(yintercept = ylim) +
  # geom_vline(xintercept = xlim) +
  geom_path(data = map, aes(x, y, group = group)) +
  ggplot2::coord_equal() +
  ggplot2::coord_equal(ylim = ylim, xlim = xlim)



ylim <- c(-4.5e6, -2e6)
xlim <- c(-0.3e6, 2.5e6)

data_sub <- data[x_m %between% xlim & y_m %between% ylim]

data_sub[, c("lon", "lat") := proj4::project(list(x_m, y_m), map_proj, inverse = TRUE)]


data_sub[, z_norm := scales::rescale(z, to = c(0, 1))]


ggplot2::ggplot() +
  
  scattermore::geom_scattermost(data_sub[!is.na(z_norm), .(lon, lat)], interpolate = TRUE, pixels = c(512, 512),
                                viridisLite::viridis(100)[1 + 99*data_sub[!is.na(z_norm), z_norm]]) +
  geom_path(data = map, aes(long, lat, group = group)) +
  ggplot2::coord_equal(ylim = c(-45, -20), xlim = c(-77, -40))

