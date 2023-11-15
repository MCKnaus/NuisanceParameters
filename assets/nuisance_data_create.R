# package data
library(DoubleML)
library(dplyr)
n = 5000
p = 10
alpha = 0.2
nuisance_data = DoubleML::make_plr_CCDDHNR2018(
  n_obs = n,
  dim_x = p,
  alpha = alpha,
  return_type = "data.frame"
)
nuisance_data$d = ifelse(nuisance_data$d > 0, 1, 0)
nuisance_data = nuisance_data %>%
  rename(w = d) %>%
  rename_with(tolower) %>%
  relocate(y, w) %>%
  as.matrix()
save(nuisance_data , file = "data/nuisance_data.RData", version = 2)
