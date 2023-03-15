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

Ks <- c(5000)
Kdf <- data.frame(matrix(Ks,nrow = 1))
names(Kdf) <- paste0("K=",Ks)
allos <- with(fhosp1 %>% filter(model == "COVIDhub-ensemble"), allocate(
      F = F,
      Q = Q,
      w = 1,
      K = 5000,
      kappa = 1,
      alpha = 1,
      dg = 1,
      eps_K = .01,
      eps_lam = .001,
      Trace = FALSE
    ))
ggplot() + xlim(0,500) +map(1:length(allos$meb), ~geom_function(fun = allos$meb[[.]]))

with(fhosp1 %>% filter(model == "CMU-TimeSeries"), oracle_allocate(
  y = value,
  w = 1,
  K = 1000,
  kappa = 1,
  alpha = 1,
  dg = 1,
  eps_K = .01,
  eps_lam = .001,
  Trace = FALSE
))

(ascores <- fhosp1 %>%
    bind_cols(Kdf) %>%
    group_by(model) %>%
  summarise(ytot = sum(value), across(starts_with("K"), ~alloscore(
  y = value,
  F = F,
  Q = Q,
  w = 1,
  K = unique(.x),
  kappa = 1,
  alpha = 1,
  dg = 1,
  eps_K = .01,
  eps_lam = .001,
  against_oracle = TRUE
))))

plot(fhosp1$F[[5]], xlim = c(0,1000))

ps8 <- fhosp1 %>% slice(8) %>% pull(ps) %>% .[[1]]
qs8 <- fhosp1 %>% slice(8) %>% pull(qs) %>% .[[1]]
cdf1 <- make_p_fn(ps = ps8, qs = qs8)

plots <- fhosp1 %>% filter(model == "CMU-TimeSeries") %>% split(.$full_location_name) %>%
  map(function(df) {
    qs <- df %>% pull(qs) %>% .[[1]]
    p <- ggplot() + geom_point(
      data = tibble(
        qs = qs,
        ps = df %>% pull(ps) %>% .[[1]]
      ),
      mapping = aes(x=qs, y=ps)
    ) + geom_line(
      data = tibble(
        x = seq(from = -1, to = 1.5*max(qs), length.out = 1000),
        y = df$F[[1]](x)
      ),
      mapping = aes(x=x,y=y)
    )
    return(p)
  })

