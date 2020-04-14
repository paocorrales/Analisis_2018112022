
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
ncfile <- nc_open("/home/paola.corrales/datosmunin/DA/DA_DATA/GOES16/OR_ABI-L1b-RadF-M3C07_G16_s20183261200367_e20183261211145_c20183261211181.nc")

proj_info <- ncatt_get(ncfile, "goes_imager_projection")
lon_0 <- proj_info$longitude_of_projection_origin
h <- proj_info$perspective_point_height + proj_info$semi_major_axis
sweep <- proj_info$sweep_angle_axis

#https://makersportal.com/blog/2018/11/25/goes-r-satellite-latitude-and-longitude-grid-projection-algorithm
# GOES-R projection info and retrieving relevant constants

#proj_info = g16nc.variables['goes_imager_projection']
#lon_origin = proj_info.longitude_of_projection_origin
#H = proj_info.perspective_point_height+proj_info.semi_major_axis

r_eq <- proj_info$semi_major_axis
r_pol <- proj_info$semi_minor_axis

# Data info
x_atr <- ncatt_get(ncfile, "x")
x <- ncvar_get(ncfile, "x")*x_atr$scale_factor + x_atr$add_offset
y_atr <- ncatt_get(ncfile, "y")
y <- ncvar_get(ncfile, "y")*y_atr$scale_factor + y_atr$add_offset


x <- c(x)
y <- c(y)

# lat/lon calc routine from satellite radian angle vector
lambda_0 <-  (lon_0*pi)/180.0

a <- sin(x)^2 + cos(x)^2*(cos(y)^2 + (sin(y)*(r_eq/r_pol))^2)
b <- -2*h*cos(x)*cos(y)
c <- h^2 - r_eq^2

r_s <- (- b - sqrt(b^2 - 4*a*c))/2*a

s_x <- r_s*cos(x)*cos(y)
s_y <- -r_s*sin(x)
s_z <- r_s*cos(x)*sin(y)

lat <- (atan((r_eq/r_pol)^2 * s_z / sqrt((h - s_x)^2 + s_y^2)))*180/pi
lon <- (lambda_0 - atan(s_y/(h - s_x)))*180/pi

