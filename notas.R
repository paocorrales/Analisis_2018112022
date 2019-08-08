obs[, ratio := obs.guess/max(1.3, min(5.6, rerr, na.rm = TRUE), na.rm = TRUE)] %>% 
  .[usage.flag == 1 & var == "t"] %>% 
  ggplot(aes(ratio)) +
  geom_density()
