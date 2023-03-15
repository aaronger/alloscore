library(tidyverse)
library(covidHubUtils)
library(distfromq)
library(alloscore)
devtools::load_all()

#hub_repo_path <- "~/research/epi/covid/covid19-forecast-hub/"
hub_repo_path = "~/covid/covid19-forecast-hub/"

inc_hosp_targets <- paste(0:30, "day ahead inc hosp")
forecasts_hosp <- load_forecasts(
#  models = c("COVIDhub-ensemble", "CMU-TimeSeries"),
#   dates = "2022-11-08",
  dates = "2021-12-27",
  date_window_size = 6,
  #locations = "US",
  types = c("quantile"),
  targets = inc_hosp_targets,
  source = "local_hub_repo",
  hub_repo_path = hub_repo_path,
  verbose = FALSE,
  as_of = NULL,
  hub = c("US")
)

forecasts_hosp <- covidHubUtils::align_forecasts(forecasts = forecasts_hosp)

truth <- load_truth(
  truth_source = "HealthData",
  target_variable = "inc hosp"
)

trad_scores <- covidHubUtils::score_forecasts(
    forecasts = forecasts_hosp %>%
        filter(
            relative_horizon == 14,
            location < 57),
    truth = truth,
    use_median_as_point = TRUE
)

fhosp1 <- forecasts_hosp %>%
  filter(
    relative_horizon == 14,
    location < 57) %>%
  select(-type) %>%
  nest(ps = quantile, qs = value) %>%
  relocate(ps, qs) %>%
  mutate(
    ps = map(ps, deframe),
    qs = map(qs, deframe)
  )

fhosp1 %>% count(model)

mkeep <- c("BPagano-RtDriven", "COVIDhub-4_week_ensemble", "COVIDhub-baseline",
    "CU-select", "IHME-CurveFit", "JHUAPL-Bucky", "JHUAPL-Gecko",
    "MUNI-ARIMA", "USC-SI_kJalpha", "UVA-Ensemble")
fhosp1 <- fhosp1 %>%
  dplyr::filter(model %in% mkeep)

fhosp1 <- fhosp1 %>%
  add_pdqr_funs(dist = "distfromq", types = c("p", "q")) %>%
  relocate(dist, F, Q)

fhosp1 <- fhosp1 %>% left_join(
  truth %>% select(location, target_end_date, value),
  by = c("location", "target_end_date"))

ys <- fhosp1 %>% filter(model == "COVIDhub-baseline") %>% select(value)

Ks <- c(seq(from = 5500, to = 30000, by = 500),seq(from = 30000, to = 120000, by = 2000))
Kdf <- data.frame(matrix(Ks,nrow = 1))
names(Kdf) <- paste0("K=",Ks)

Ks_large <- seq(from = 30000, to = 120000, by = 2000)
K_large_df <- data.frame(matrix(Ks_large,nrow = 1))
names(K_large_df) <- paste0("K=",Ks_large)

(ascores <- fhosp1 %>%
    bind_cols(Kdf) %>%
    group_by(model) %>%
  summarise(
    ytot = sum(value),
    across(starts_with("K"), ~alloscore:::alloscore(
        y = value,
        F = F,
        Q = Q,
        w = 1,
        K = unique(.x),
        kappa = 1,
        alpha = 1,
        dg = 1,
        eps_K = .01,
        eps_lam = 1e-3,
        against_oracle = TRUE,
        verbose = TRUE
        ))
  ))


ascores_long <- cbind(ascores) %>%
    tidyr::pivot_longer(
        starts_with("K="),
        names_prefix = "K=",
        names_to = "K"
    ) %>%
    dplyr::mutate(K = as.numeric(K))

p <- ggplot(
    data = ascores_long,
    mapping = aes(x = as.numeric(K), y = value, color = model, linetype = model)
) +
    geom_line() +
    theme_bw()

library(plotly)
ggplotly(p)

ao <- alloscore::oracle_allocate(ys$value,
                                 w = 1, K = 10000, alpha = 1,
                                 eps_K = .001, eps_lam = 1e-4, Trace = TRUE)

trad_scores %>%
    filter(model %in% mkeep) %>%
    group_by(model) %>%
    summarize(wis = sum(wis)) %>%
    arrange(wis)


for (m in unique(fhosp1$model)) {
    m <- "USC-SI_kJalpha"
    m <- "MUNI-ARIMA"
    print(m)
    alloscore:::alloscore(
        y = fhosp1 %>% filter(model == m) %>% pull(value),
        F = fhosp1 %>% filter(model == m) %>% pull(F),
        Q = fhosp1 %>% filter(model == m) %>% pull(Q),
        w = 1,
        K = 38555,
        kappa = 1,
        alpha = 1,
        dg = 1,
        eps_K = .001,
        eps_lam = 1e-5,
        against_oracle = FALSE
    )
    alloscore:::alloscore(
        y = fhosp1 %>% filter(model == m) %>% pull(value),
        F = fhosp1 %>% filter(model == m) %>% pull(F),
        Q = fhosp1 %>% filter(model == m) %>% pull(Q),
        w = 1,
        K = 38554,
        kappa = 1,
        alpha = 1,
        dg = 1,
        eps_K = .001,
        eps_lam = 1e-5,
        against_oracle = FALSE
    )
    a <- alloscore:::allocate(
        F = fhosp1 %>% filter(model == m) %>% pull(F),
        Q = fhosp1 %>% filter(model == m) %>% pull(Q),
        w = 1,
        K = 38554,
        kappa = 1,
        alpha = 1,
        dg = 1,
        eps_K = .001,
        eps_lam = 1e-5,
        Trace = TRUE
    )
    a2 <- alloscore:::allocate(
      F = fhosp1 %>% filter(model == m) %>% pull(F),
      Q = fhosp1 %>% filter(model == m) %>% pull(Q),
      w = 1,
      K = 38555,
      kappa = 1,
      alpha = 1,
      dg = 1,
      eps_K = .001,
      eps_lam = 1e-5,
      Trace = TRUE
    )

    # fhosp1 %>% filter(model == m) %>% slice(1) %>% pull(ps)
    # fhosp1 %>% filter(model == m) %>% slice(1) %>% pull(qs)
    # Q <- fhosp1 %>% filter(model == m) %>% slice(1) %>% pull(Q) %>% `[[`(1)

    # fhosp1 %>% filter(model == m) %>% as.data.frame() %>% nrow()
}

conv_dat <- a$xs %>% map(enframe) %>% enframe() %>%
  rename(tau = name) %>% unnest(cols = value) %>% mutate(location = as.factor(name)) %>% mutate(K="38554")
conv_dat2 <- a2$xs %>% map(enframe) %>% enframe() %>%
  rename(tau = name) %>% unnest(cols = value) %>% mutate(location = as.factor(name))%>% mutate(K="38555")
df <- bind_rows(conv_dat, conv_dat2)

ggplot(data=df %>% filter(tau>1), aes(x=tau,y=value,color=location,linetype=as.factor(K)))+geom_line()

ggplot() + xlim(c(12, 15)) + geom_function(fun = function(x) pnorm(x, log.p = TRUE), n=1000)

  map(1:51, ~geom_function(fun = a$meb[[.]]))


plot_by_tau <- function(allo_trace) {
  dat <- map(allo_trace$xis,~enframe(.)) %>% enframe() %>% unnest()
}



plot(fhosp1$F[[8]], xlim = c(0,500))

ps8 <- fhosp1 %>% slice(8) %>% pull(ps) %>% .[[1]]
qs8 <- fhosp1 %>% slice(8) %>% pull(qs) %>% .[[1]]
cdf1 <- make_p_fn(ps = ps8, qs = qs8)

plots <- fhosp1 %>% filter(model == "CMU-TimeSeries") %>% split(.$full_location_name) %>%
  map(function(df) {
    qs <- df %>% pull(qs) %>% .[[1]]
    p <- ggplot(
      data.frame(
        x=seq(from = -1, to = 1.5*max(qs), length.out = 1000), y = ))
  })

i <- 2
i <- i + 1
ggplot() + geom_function(fun=fhosp1$F[[i]], xlim = c(0, 1000))
