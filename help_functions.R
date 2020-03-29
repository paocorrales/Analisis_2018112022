# Help functions


# read_diag ---------------------------------------------------------------
# Read diagfiles and tidy uv observations

read_diag <- function(path, variable) {
  library(foreach)
  library(doParallel)
  
  files <- Sys.glob(paste0(path, ".mem*"))
  
  myCluster <- makeCluster(20)
  registerDoParallel(myCluster)
  
  obs <- foreach(f = 1:length(files),
                 .packages = c("data.table", "metR", "lubridate", "tidyverse"),
                 .export = c("files", "variable"),
                 .combine = "rbind") %dopar% {
                   
                   # sink("log.txt", append=TRUE)
                   
                   diag <- fread(files[f]) %>% 
                     .[V10 == 1] %>% 
                     .[, exp := basename(dirname(files[f]))] %>% 
                     .[, mem := str_extract(files[f], "\\d{3}$")] %>% 
                     .[, date := ymd_hms(str_extract(files[f], "\\d{14}"))] %>% 
                     .[, c("V2", "V4") := NULL]
                   
                   # cat("Archivo ", basename(files[f]))
                   
                   colnames(diag) <- c("var", "stationID", "type", "dhr", "lat", "lon", "pressure", "usage.flag", "flag.prep", "obs", "obs.guess", "obs2", "obs.guess2", "rerr", "exp", "mem", "date")
                   
                   if ("uv" %in% variable & length(variable) == 1) {
                     
                     diag <- diag[var == "uv"] %>% 
                       melt(measure.vars = c("obs", "obs2", "obs.guess", "obs.guess2")) %>% 
                       .[, var := if_else(str_detect(variable, "2"), "v", "u")] %>% 
                       .[, variable := str_remove(variable, "2")] 
                     
                     vars <- rlang::syms(setdiff(names(diag), "value")) 
                     diag <- diag %>% 
                       distinct(!!!vars, .keep_all = TRUE) %>%  
                       pivot_wider(names_from = variable, values_from = value) %>% 
                       setDT %>% 
                       .[, id := 1:.N, by = mem] 
                     
                   } else if ("uv" %in% variable & length(variable) != 1) {
                     
                     
                     uv <- diag[var == "uv"] %>% 
                       melt(measure.vars = c("obs", "obs2", "obs.guess", "obs.guess2")) %>% 
                       .[, var := if_else(str_detect(variable, "2"), "v", "u")] %>% 
                       .[, variable := str_remove(variable, "2")] 
                     
                     vars <- rlang::syms(setdiff(names(uv), "value")) 
                     uv <- uv %>% 
                       distinct(!!!vars, .keep_all = TRUE) %>% 
                       pivot_wider(names_from = variable, values_from = value) %>% 
                       setDT
                     
                     variable <- c(variable, "u", "v")
                     
                     diag <- diag[var != "uv", -c("obs2", "obs.guess2"), with = FALSE] %>% 
                       rbind(uv) %>% 
                       .[var %in% variable] %>% 
                       .[, id := 1:.N, by = mem]  
                     
                   } else {
                     diag <- diag[var %in% variable, -c("obs2", "obs.guess2"), with = FALSE] %>% 
                       .[, id := 1:.N, by = mem] 
                   }
                   
                   diag[, obs := ifelse(obs == -1e+05, NA, obs)][]
                 }
  stopCluster(myCluster)
  return(obs)
}

# input_obs_error ---------------------------------------------------------

input_obs_error <- function(variable, type, nivel, path_to_errtable = "errtable.csv") {
  errtable <- fread(path_to_errtable) %>% 
    .[, ":="(u = uv,
             v = uv, 
             uv = NULL,
             pw = NULL)] %>% 
    melt(id.vars = c("type", "nivel"))
  
  
  data.table(variable = variable, type = type, pressure = nivel) %>% 
    .[, nivel := errtable[metR::Similar(nivel, pressure), unique(nivel)], by = pressure] %>% 
    errtable[., on = .NATURAL] %>% 
    .[, .(value, nivel)]
  
}


# CR_bias -----------------------------------------------------------------

get_cr <-  function(dt, tipo = "superficie") {
  
  if ("temporal" %in% tipo) {
    if ("superficie" %in% tipo) {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, guess := obs - obs.guess] %>% 
        .[, ":="(guess.mean = mean(guess),
                 obs.guess.mean = mean(obs.guess),                # media de la innovación (y0 - H(xf)) para el ensamble 
                 var.guess = var(guess)), by = .(id, date)] %>%   # varianza de H(xf) calculada sobre el ensamble
        .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type, date)] %>% # innovación promediada sobre el dominio para cada tipo de obs
        .[, .(rmsi = mean((obs.guess.mean - d.mean)^2),                       # <(d - <d>)^2>
              spread = error[1]^2 + mean(var.guess, na.rm = TRUE),            # error^2 + media de la varianza de H(xf) sobre el dominio
              bias = mean(obs.guess, na.rm = TRUE),                           # bias hecho y derecho, obs.guess = y0 - H(xf)
              count = .N),  by = .(var, type, date)] %>%                      # Cantidad de observaciones por tipo (x60)
        .[, cr := spread / rmsi]
    } else {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, guess := obs - obs.guess] %>% 
        .[, ":="(guess.mean = mean(guess),
                 obs.guess.mean = mean(obs.guess),
                 var.guess = var(guess)), by = .(id, date)] %>% 
        .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type, nivel.error, date)] %>% 
        .[, .(rmsi = mean((obs.guess.mean - d.mean)^2),
              spread = error[1]^2 + mean(var.guess, na.rm = TRUE),
              bias = mean(obs.guess, na.rm = TRUE),
              count = .N),  by = .(var, type, nivel.error, date)] %>% 
        .[, cr := spread / rmsi]
    } 
  } else {
    if ("superficie" %in% tipo) {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, guess := obs - obs.guess] %>% 
        .[, ":="(guess.mean = mean(guess),
                 obs.guess.mean = mean(obs.guess),
                 var.guess = var(guess)), by = .(id, date)] %>% 
        .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type)] %>% 
        .[, .(rmsi = mean((obs.guess.mean - d.mean)^2),
              spread = error[1]^2 + mean(var.guess, na.rm = TRUE),
              bias = mean(obs.guess, na.rm = TRUE),
              count = .N),  by = .(var, type)] %>% 
        .[, cr := spread / rmsi]
    } else {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, guess := obs - obs.guess] %>% 
        .[, ":="(guess.mean = mean(guess),
                 obs.guess.mean = mean(obs.guess),
                 var.guess = var(guess)), by = .(id, date)] %>% 
        .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type, nivel.error)] %>% 
        .[, .(rmsi = mean((obs.guess.mean - d.mean)^2),
              spread = error[1]^2 + mean(var.guess, na.rm = TRUE),
              bias = mean(obs.guess, na.rm = TRUE),
              count = .N),  by = .(var, type, nivel.error)] %>% 
        .[, cr := spread / rmsi]
    }
  }
}


# RCRV --------------------------------------------------------------------

get_RCRV <- function(dt, tipo = "superficie") {
  
  if ("temporal" %in% tipo) {
    if ("superficie" %in% tipo) {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, ":="(mean.guess = mean(obs - obs.guess, na.rm = TRUE),
                 sd.guess = sd(obs - obs.guess, na.rm = TRUE)), by = .(id, date)] %>% 
        .[, y := (obs - mean.guess)/sqrt(sd.guess^2 + error^2), by = .(id, date)] %>% 
        .[, .(mean.y = mean(y, na.rm = TRUE),
              sd.y = sd(y, na.rm = TRUE)), by = .(var, type, date)]
    } else {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, ":="(mean.guess = mean(obs - obs.guess, na.rm = TRUE),
                 sd.guess = sd(obs - obs.guess, na.rm = TRUE)), by = .(id, date)] %>% 
        .[, y := (obs - mean.guess)/sqrt(sd.guess^2 + error^2), by = .(id, date)] %>% 
        .[, .(mean.y = mean(y, na.rm = TRUE),
              sd.y = sd(y, na.rm = TRUE)), by = .(var, type, nivel.error, date)]
    } 
  } else {
    if ("superficie" %in% tipo) {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, ":="(mean.guess = mean(obs - obs.guess, na.rm = TRUE),
                 sd.guess = sd(obs - obs.guess, na.rm = TRUE)), by = .(id, date)] %>% 
        .[, y := (obs - mean.guess)/sqrt(sd.guess^2 + error^2), by = .(id, date)] %>% 
        .[, .(mean.y = mean(y, na.rm = TRUE),
              sd.y = sd(y, na.rm = TRUE)), by = .(var, type)]
    } else {
      dt[usage.flag == 1 & !is.na(obs)] %>% 
        .[, ":="(mean.guess = mean(obs - obs.guess, na.rm = TRUE),
                 sd.guess = sd(obs - obs.guess, na.rm = TRUE)), by = .(id, date)] %>% 
        .[, y := (obs - mean.guess)/sqrt(sd.guess^2 + error^2), by = .(id, date)] %>% 
        .[, .(mean.y = mean(y, na.rm = TRUE),
              sd.y = sd(y, na.rm = TRUE)), by = .(var, type, nivel.error)]
    }
  }
  
}


# Label box ---------------------------------------------------------------

cut_round <- function(x, breaks) {
  labels <- na.omit((breaks + data.table::shift(breaks, -1))/2)
  cuts <- cut(x, breaks = breaks, labels = labels)
  
  as.numeric(as.character(cuts)) 
}


# Wrap FSS ----------------------------------------------------------------

# Usa la función fss del paquete verification pero previamente reorganiza las
# variables en matrices. También puede iterar para distintos q (valor de pp) y
# w (tamaño de la caja = w2+1)

FSS <- function(lon, lat, obs, fcst, q, w){
  dt <- data.table(obs.binary = obs,
                   fcst.binary = fcst, 
                   lon,
                   lat)
  
  fcst <- dt %>%
    dcast(lon ~ lat, value.var = "fcst.binary") %>%
    .[, -1] %>% 
    as.matrix()
  
  obs <- dt %>%
    dcast(lon ~ lat, value.var = "obs.binary") %>%
    .[, -1] %>% 
    as.matrix()
  
  out <- purrr::map_dfr(q, function(q) {
    # browser()
    fcst_q <- fcst >= q
    obs_q <- obs >= q
    
    return <- list(fss = purrr::map_dbl(w, ~ verification::fss(obs_q, fcst_q, .x)),
         w = w,
         q = rep(q, length(w)))
    message(paste("Listo q = ", q))
    return(return)
  })
  
  return(out)
}


read_radiosonde_relampago <- function(file){
  # Leo línea por línea
  lines <- readLines(file)
  
  # Indices donde comienza cada sondeo
  idx <- which(grepl("Data Type:", lines))
  idx <- c(idx, length(temp)+1)
  
  soundings <- list()
  for (i in seq_len(length(idx)-1)) { 
    
    out <- read.table(text = lines[(idx[i] + 15):(idx[i + 1] - 1)]) %>% 
      as.data.table()
    
    names <- strsplit(lines[idx[i] + 12], " ")[[1]]
    names <- names[names != ""]
    colnames(out) <- names
    
    launch <- lubridate::ymd_hms(strsplit(lines[idx[i] + 4], "    ")[[1]][2])
    nominal_launch <- lubridate::ymd_hms(strsplit(lines[idx[i] + 11], "):")[[1]][2])
    site <- strsplit(lines[idx[i] + 2], "         ")[[1]][2]  
    
    out[, ":="(Site = site,
               nominal_launch_time = nominal_launch,
               launch_time = launch)]
    
  }
}