# Calculate vertical profile of T and Q

library(tidyverse)
library(data.table)
library(metR)
library(lubridate)


era_path <- "/home/paola.corrales/datosmunin3/DATA_EXP/"
exp_path <- "/home/paola.corrales/datosmunin3/EXP/"

lat_band <- c(40, 299)
lon_band <- c(70, 180)

exp <- "E6_long"
ini_date <- ymd_hms("20181108000000")
ciclos <- 379

dates <- seq.POSIXt(ini_date, by = "hour",
                    length.out = ciclos)

levs <- c(seq(1000, 500, length.out = 20), seq(500, 50, by = -50))
levs <- unique(levs)

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

rh <- function(QVAPOR, P, T) {
  P <- P*0.01      # Debe estar en hPa
  T <- T - 273.15  # Debe estar en Celsius
  es <- 6.112 * exp(17.67*(T)/(T + 273.15 - 29.65))

  qvs <- es/(P - (1 - 0.622)*es)

  rh <- 100*pmax(pmin(QVAPOR/qvs, 1), 0)

  return(rh)
}

for (d in 1:length(dates)) {
  print(paste("met_em", dates[d]))

  gue_file <- paste0(exp_path, exp, "/GUESS/", format(dates[d], "%Y%m%d%H%M%S"), "/wrfarw.ensmean")
  
  if (!file.exists(gue_file)) {
    message(paste0("Archivo para ", dates[d], " no existe"))
    next
  }

  gue <- ReadNetCDF(gue_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, RH := rh(QVAPOR, T, P)] %>%
    .[, c("U", "V") := uvmet(gue_file)] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>%
    .[, date := dates[d]] %>%
    .[south_north %between% lat_band & west_east %between% lon_band] %>%
    .[, .(lev = levs,
          T = approx(P, T, xout = levs*100)$y,
          RH = approx(P, RH, xout = levs*100)$y,
          QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
          U = approx(P, U, xout = levs*100)$y,
          V = approx(P, V, xout = levs*100)$y,
          SPD = approx(P, SPD, xout = levs*100)$y),
      by = .(date, south_north, west_east)]

  ana_file <- paste0(exp_path, exp, "/ANA/", format(dates[d], "%Y%m%d%H%M%S"), "/analysis.ensmean")

  ana <- ReadNetCDF(ana_file, vars = c("P", "PB", "T", "QVAPOR")) %>%
    .[, P := P + PB] %>%
    .[, T := (T + 300)*(P/100000)^(2/7)] %>%
    .[, RH := rh(QVAPOR, T, P)] %>%
    .[, c("U", "V") := uvmet(ana_file)] %>%
    .[, SPD := sqrt(U^2 + V^2)] %>%
    .[, date := dates[d]] %>%
    .[south_north %between% lat_band & west_east %between% lon_band] %>%
    .[, .(lev = levs,
          T = approx(P, T, xout = levs*100)$y,
          RH = approx(P, RH, xout = levs*100)$y,
          QVAPOR = approx(P, QVAPOR, xout = levs*100)$y,
          U = approx(P, U, xout = levs*100)$y,
          V = approx(P, V, xout = levs*100)$y,
          SPD = approx(P, SPD, xout = levs*100)$y),
      by = .(date, south_north, west_east)]

  # Medias

  # era_media <- era %>%
  #    .[, .(T = mean(T, na.rm = TRUE),
  #          RH = mean(RH, na.rm = TRUE),
  #          QVAPOR = mean(QVAPOR, na.rm = TRUE),
  #          U = mean(U, na.rm = TRUE),
  #          V = mean(V, na.rm = TRUE),
  #          SPD = mean(SPD, na.rm = TRUE)), by = .(lev, date)]
  #
  # ana_media <- ana %>%
  # .[, .(T = mean(T, na.rm = TRUE),
  #       RH = mean(RH, na.rm = TRUE),
  #       QVAPOR = mean(QVAPOR, na.rm = TRUE),
  #       U = mean(U, na.rm = TRUE),
  #       V = mean(V, na.rm = TRUE),
  #       SPD = mean(SPD, na.rm = TRUE)), by = .(lev, date)]

  # RMSE

  out <- gue %>%
    melt(id = c("date", "south_north", "west_east", "lev"), value.name = "gue") %>%
    .[ana %>% melt(id = c("date", "south_north", "west_east", "lev"), value.name = "ana"), on = .NATURAL] %>%
    .[, .(rmse = sd(ana - gue, na.rm = TRUE),
          bias = mean(ana - gue, na.rm = TRUE),
          ana = mean(ana, na.rm = TRUE),
          gue = mean(gue, na.rm = TRUE)), by = .(date, lev, variable)] %>%
    .[, exp := exp]

  if (d == 1) {
    # out_guess <- guess
    out_ana <- out
  } else {
    # out_guess <- rbind(out_guess, guess)
    out_ana <- rbind(out_ana, out)
  }

}

# fwrite(out_guess, paste0(path, "/derived_data/perfiles_guess_", exp, ".csv"))
fwrite(out_ana, paste0("~/datosmunin3/EXP/derived_data/fcst_det/perfiles_ana_gue_", exp, ".csv"))

