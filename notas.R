obs[, ratio := obs.guess/max(1.3, min(5.6, rerr, na.rm = TRUE), na.rm = TRUE)] %>% 
  .[usage.flag == 1 & var == "t"] %>% 
  ggplot(aes(ratio)) +
  geom_density()


diag <- fread("E1/test/test.csv", na.strings = c("0.100E+11", "-0.100E+06", "-99999.90", "-100000.00"))
