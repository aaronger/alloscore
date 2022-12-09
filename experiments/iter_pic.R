library(tidyverse)

(dat <- tibble(
  alpha = c(.2,.6,.3,.9),
  w = c(3,2,3,5),
  dists_and_params = list(
    list(dist = "norm", mean = 3, sd = 2),
    list(dist = "norm", mean = 7, sd = 9),
    list(dist = "exp", rate = .3),
    list(dist = "gamma", shape = 1, scale = 4)
  ),
  name = map(dists_and_params,
             ~toString(paste0(.[["dist"]], ", ",
                              toString(paste(names(.)[-1], .[-1], sep = "=")))))) %>%
  unnest_wider(col = dists_and_params) %>%
  add_pdqr_funs()
)

(dat_allo <- dat %>%
    mutate(
      q_alpha = map2_dbl(Q, alpha, exec),
      allo = allocate(
        F = F,
        Q = Q,
        alpha = alpha,
        w = w,
        K = 10,
        eps_K = .001,
        eps_lam = .001
      )
    ) %>%
    relocate(allo)
  )

  dat_allo <- dat_allo %>% mutate(
    Lambda = pmap(., margexb_fun))

with(dat_allo, ggplot() + xlim(-1, 2) +
  map(1:4, ~ geom_function(aes(color = name[[.]]), fun = Lambda[[.]])) +
    geom_hline(yintercept = 0) +
    labs(color = "Distribution", y = expression(lambda(x)), x="x")
  )

