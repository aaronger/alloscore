library(covidHubUtils)

inc_hosp_targets <- paste(0:130, "day ahead inc hosp")
forecasts_hosp <- load_forecasts(
  models = c("COVIDhub-ensemble", "CMU-TimeSeries"),
  dates = "2022-11-08",
  date_window_size = 6,
  #locations = "US",
  types = c("quantile"),
  targets = inc_hosp_targets,
  source = "zoltar",
  verbose = FALSE,
  as_of = NULL,
  hub = c("US")
)

forecasts_hosp
