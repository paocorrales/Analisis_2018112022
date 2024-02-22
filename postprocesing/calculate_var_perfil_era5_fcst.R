# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)


era_path <- "/home/paola.corrales/datosmunin3/DATA_EXP/"
exp_path <- "/home/paola.corrales/datosmunin3/EXP/"

lat_band <- c(40, 299)
lon_band <- c(70, 180)

# lon_band = c(-66.97, -54.51)
# lat_band = c(-38.14, -20.10)

all_exp <- c("E2", "E5", "E6", "E9")
run <- "fcst"
ini_date <- ymd_hms("20181122060000")
ciclos <- 30

# dates <- seq.POSIXt(ini_date, by = "hour",
#                     length.out = ciclos)
lead_time <- seq(0, ciclos)

# levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
# levs <- unique(levs)

uvmet_era <- function(ncfile) {

  nc <- ncdf4::nc_open(ncfile)

  # https://github.com/NCAR/wrf-python/blob/b40d1d6e2d4aea3dd2dda03aae18e268b1e9291e/src/wrf/g_uvmet.py#L179
  true_lat1 <- ncdf4::ncatt_get(nc, 0, "TRUELAT1")[[2]]
  true_lat2 <- ncdf4::ncatt_get(nc, 0, "TRUELAT2")[[2]]
  cen_lon <- ncdf4::ncatt_get(nc, 0, "STAND_LON")[[2]]

  lon <- ncdf4::ncvar_get(nc, "XLONG_M")
  lat <- ncdf4::ncvar_get(nc, "XLAT_M")

  u <- ncdf4::ncvar_get(nc, "UU")
  v <- ncdf4::ncvar_get(nc, "VV")

  rpd <- pi/180.0

  if ((abs(true_lat1 - true_lat2) > 0.1) & (abs(true_lat2 - 90.) > 0.1)) {
    cone = (log(cos(true_lat1*rpd)) -
              log(cos(true_lat2*rpd)))
    cone = (cone /
              (log(tan((45.- abs(true_lat1/2.))*rpd))
               - log(tan((45.- abs(true_lat2/2.)) *
                           rpd))))
  } else {
    cone = sin(abs(true_lat1)*rpd)
  }

  # https://github.com/NCAR/wrf-python/blob/d9585354c0e2a75a0f7c1d6b200d353f5e4eb084/fortran/wrf_user.f90#L801

  longca <- lon - cen_lon

  longca[longca > 180] <- longca[longca > 180] - 360
  longca[longca < -180] <- longca[longca < -180] + 360

  longcb <- longca*cone*rpd            # Hemisferio norte
  longcb[lat < 0] <- -longcb[lat < 0]  # Corrección si el dominio está en el HS

  longca <- cos(longcb)
  longcb <- sin(longcb)

  # destagger
  u <- 0.5*(u[1:(nrow(u)-1), ,] + u[c(2:nrow(u)), , ])
  v <- 0.5*(v[, 1:(ncol(v)-1), ] + v[, c(2:ncol(v)), ])
  longcb <- rray::rray_broadcast(longcb, dim(v))
  longca <- rray::rray_broadcast(longca, dim(u))
  # desrotación
  umet <- v * longcb + u * longca
  vmet <- v * longca - u * longcb

  list(c(umet),
       c(vmet))
}

uvmet <- function(ncfile) {

  nc <- ncdf4::nc_open(ncfile)

  # https://github.com/NCAR/wrf-python/blob/b40d1d6e2d4aea3dd2dda03aae18e268b1e9291e/src/wrf/g_uvmet.py#L179
  true_lat1 <- ncdf4::ncatt_get(nc, 0, "TRUELAT1")[[2]]
  true_lat2 <- ncdf4::ncatt_get(nc, 0, "TRUELAT2")[[2]]
  cen_lon <- ncdf4::ncatt_get(nc, 0, "STAND_LON")[[2]]

  lon <- ncdf4::ncvar_get(nc, "XLONG")
  lat <- ncdf4::ncvar_get(nc, "XLAT")

  u <- ncdf4::ncvar_get(nc, "U")
  v <- ncdf4::ncvar_get(nc, "V")

  rpd <- pi/180.0

  if ((abs(true_lat1 - true_lat2) > 0.1) & (abs(true_lat2 - 90.) > 0.1)) {
    cone = (log(cos(true_lat1*rpd)) -
              log(cos(true_lat2*rpd)))
    cone = (cone /
              (log(tan((45.- abs(true_lat1/2.))*rpd))
               - log(tan((45.- abs(true_lat2/2.)) *
                           rpd))))
  } else {
    cone = sin(abs(true_lat1)*rpd)
  }

  # https://github.com/NCAR/wrf-python/blob/d9585354c0e2a75a0f7c1d6b200d353f5e4eb084/fortran/wrf_user.f90#L801

  longca <- lon - cen_lon

  longca[longca > 180] <- longca[longca > 180] - 360
  longca[longca < -180] <- longca[longca < -180] + 360

  longcb <- longca*cone*rpd            # Hemisferio norte
  longcb[lat < 0] <- -longcb[lat < 0]  # Corrección si el dominio está en el HS

  longca <- cos(longcb)
  longcb <- sin(longcb)

  # destagger
  u <- 0.5*(u[1:(nrow(u)-1), ,] + u[c(2:nrow(u)), , ])
  v <- 0.5*(v[, 1:(ncol(v)-1), ] + v[, c(2:ncol(v)), ])
  longcb <- rray::rray_broadcast(longcb, dim(v))
  longca <- rray::rray_broadcast(longca, dim(u))
  # desrotación
  umet <- v * longcb + u * longca
  vmet <- v * longca - u * longcb

  list(c(umet),
       c(vmet))
}

qv <- function(RH, T, P) {

  # qv <- RH*exp(17.67*(T - 273.16)/(T - 29.65))/0.263*P
  es <- ClausiusClapeyron(T)
  e <- (RH/100) * es
  w <- MixingRatio(P, e, epsilon = 0.622)
  qv <- w/(w+1)
  return(qv)
}

# rh <- function(QVAPOR, P, T) {
#   P <- P*0.01      # Debe estar en hPa
#   T <- T - 273.15  # Debe estar en Celsius
#   es <- 6.112 * exp(17.67*(T)/(T + 273.15 - 29.65))
#
#   qvs <- es/(P - (1 - 0.622)*es)
#
#   rh <- 100*pmax(pmin(QVAPOR/qvs, 1), 0)
#
#   return(rh)
# }

rh <- function(QVAPOR, T, P) {

  # T en kelvin, P en Pa
  # QVAPOR dimensionless

  rh <- (0.263 * P * QVAPOR) / (exp((17.67*(T - 273.16))/(T - 29.65)))
  return(rh)
}

purrr::map(all_exp, function(exp){
  print(paste("arranca ", exp))

perfiles <- map_df(lead_time, function(d) {
  print(paste("lead_time", d))

  file_wrf <- paste0(exp_path, exp, "/", toupper(run), "/",
                     format(ini_date, "%Y%m%d%H"), "/NPP/NPP_",
                     format(ini_date, "%Y"), "-",
                     format(ini_date, "%m"), "-",
                     format(ini_date, "%d"), "_",
                     format(ini_date, "%H"), "_FC",
                     formatC(d, width = 2, flag = "0"),
                     ".nc")


  if (!file.exists(file_wrf)) {
    return(NULL)
  }

  fcst <- ReadNetCDF(file_wrf, vars = c("T", QVAPOR = "Q", "U", "V")) %>%
    .[lat %between% lat_band & lon %between% lon_band] %>%
    .[, QVAPOR := QVAPOR/1000] %>%
    .[, RH := rh(QVAPOR, T, lev*100)] %>%
    .[, ":="(exp = exp,
             date =  ini_date + hours(d))] %>%
    .[, .(T = mean(T, na.rm = TRUE),
          RH = mean(RH, na.rm = TRUE),
          QVAPOR = mean(QVAPOR, na.rm = TRUE),
          U = mean(U, na.rm = TRUE),
          V = mean(V, na.rm = TRUE)), by = .(lev, date, lon, lat)]

  levs <- unique(fcst$lev)

  era_file <- paste0(era_path, "ERA5/met_em.d01.", format(ini_date  + hours(d), "%Y-%m-%d_%H:%M:%S"), ".nc")

  era <- ReadNetCDF(era_file, vars = c(P = "PRES", T = "TT", "RH")) %>%
    .[, QVAPOR := qv(RH, T, P)] %>%
    .[, P := P/100] %>%
    .[, c("U", "V") := uvmet_era(era_file)] %>%
    .[, date := ini_date + hours(d)] %>%
    .[south_north %between% lat_band & west_east %between% lon_band] %>%
    .[, .(lev = levs,
          T = approx(P, T, xout = levs)$y,
          RH = approx(P, RH, xout = levs)$y,
          QVAPOR = approx(P, QVAPOR, xout = levs)$y,
          U = approx(P, U, xout = levs)$y,
          V = approx(P, V, xout = levs)$y),
      by = .(date, lat = south_north, lon = west_east)]



  out <- era %>%
    melt(id = c("date", "lat", "lon", "lev"), value.name = "era") %>%
    .[fcst %>% melt(id = c("date", "lat", "lon", "lev"), value.name = "fcst"), on = .NATURAL] %>%
    .[, .(rmse = sd(fcst - era, na.rm = TRUE),
          bias = mean(fcst - era, na.rm = TRUE),
          fcst = mean(fcst, na.rm = TRUE),
          era = mean(era, na.rm = TRUE)), by = .(date, lev, variable)] %>%
    .[, exp := exp]


})

# fwrite(out_guess, paste0(path, "/derived_data/perfiles_guess_", exp, ".csv"))
fwrite(perfiles, paste0("~/datosmunin3/EXP/derived_data/perfiles_fcst_", exp, "_", format(ini_date, "%Y%m%d%H"), ".csv"))

exp
})
