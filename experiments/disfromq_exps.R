library(tidyverse)
library(covidHubUtils)
library(distfromq)

inc_hosp_targets <- paste(0:130, "day ahead inc hosp")
forecasts_hosp <- load_forecasts(
  models = c("COVIDhub-ensemble", "CMU-TimeSeries"),
  dates = "2022-11-08",
  date_window_size = 6,
  #locations = "US",
  types = c("quantile"),
  targets = inc_hosp_targets,
  source = "local_hub_repo",
  hub_repo_path = "~/covid/covid19-forecast-hub/",
  verbose = FALSE,
  as_of = NULL,
  hub = c("US")
)

fhosp1 <- forecasts_hosp %>% filter(
  #model == "COVIDhub-ensemble",
  horizon == 14,
  location < 57) %>% select(-type) %>%
  nest(ps = quantile, qs = value) %>%
  relocate(ps, qs) %>%
  mutate(
    ps = map(ps, deframe),
    qs = map(qs, deframe)
    )

fhosp1 <- fhosp1 %>% add_pdqr_funs(dist = "distfromq", types = c("p", "q")) %>%
  relocate(dist, F, Q)

truth <- load_truth(
  truth_source = "HealthData",
  target_variable = "inc hosp"
)

fhosp1 <- fhosp1 %>% left_join(
  truth %>% select(location, target_end_date, value),
  by = c("location", "target_end_date"))

Ks <- c(1000)
Kdf <- data.frame(matrix(Ks,nrow = 1))
names(Kdf) <- paste0("K=",Ks)

(ascores <- fhosp1 %>%
    bind_cols(Kdf) %>%
    group_by(model) %>%
  summarise(ytot = sum(value), across(starts_with("K"), ~alloscore(
  y = value,
  F = F,
  Q= Q,
  w = 1,
  K = unique(.x),
  kappa = 1,
  alpha = 1,
  dg = 1,
  eps_K = .01,
  eps_lam = .01,
  against_oracle = FALSE
))))

plot(fhosp1$F[[8]], xlim = c(40,55))

