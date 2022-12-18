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
  model == "COVIDhub-ensemble",
  horizon == 7,
  location < 57) %>% select(-type) %>%
  nest(ps = quantile, qs = value) %>%
  relocate(ps, qs) %>%
  mutate(
    ps = map(ps, deframe),
    qs = map(qs, deframe)
    )

fhosp1 <- fhosp1 %>% add_pdqr_funs(dist = "distfromq")




