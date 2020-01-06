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
        .[, variable := str_remove(variable, "2")] %>% 
        pivot_wider(names_from = variable, values_from = value) %>% 
        setDT %>% 
        .[, id := 1:.N, by = mem] 
    
  } else if ("uv" %in% variable & length(variable) != 1) {
    uv <- diag[var == "uv"] %>% 
      melt(measure.vars = c("obs", "obs2", "obs.guess", "obs.guess2")) %>% 
      .[, var := if_else(str_detect(variable, "2"), "v", "u")] %>% 
      .[, variable := str_remove(variable, "2")] %>% 
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
  
  if (tipo == "superficie") {
    dt[usage.flag == 1 & !is.na(obs)] %>% 
      .[, guess := obs - obs.guess] %>% 
      .[, ":="(guess.mean = mean(guess),
               obs.guess.mean = mean(obs.guess),
               var.guess = var(guess)), by = .(id, date)] %>% 
      .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type)] %>% 
      .[, .(dd = mean((obs.guess.mean - d.mean)^2),
            rmse = error[1]^2 + mean(var.guess, na.rm = TRUE),
            bias = mean(obs.guess, na.rm = TRUE),
            count = .N),  by = .(var, type)] %>% 
      .[, cr := rmse / dd]
  } else {
    dt[usage.flag == 1 & !is.na(obs)] %>% 
      .[, guess := obs - obs.guess] %>% 
      .[, ":="(guess.mean = mean(guess),
               obs.guess.mean = mean(obs.guess),
               var.guess = var(guess)), by = .(id, date)] %>% 
      .[, ":="(d.mean = mean(obs.guess.mean)), by = .(var, type, nivel.error)] %>% 
      .[, .(dd = mean((obs.guess.mean - d.mean)^2),
            rmse = error[1]^2 + mean(var.guess, na.rm = TRUE),
            bias = mean(obs.guess, na.rm = TRUE),
            count = .N),  by = .(var, type, nivel.error)] %>% 
      .[, cr := rmse / dd]
  }
}
