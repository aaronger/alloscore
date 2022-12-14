library(tidyverse)
library(patchwork)
library(alloscore)

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
K=10
(dat_allo <- dat %>%
    mutate(
      q_alpha = map2_dbl(Q, alpha, exec),
      allo = allocate(
        F = F,
        Q = Q,
        alpha = alpha,
        w = w,
        K = K,
        eps_K = .001,
        eps_lam = .001
      )
    ) %>%
    relocate(allo)
  )

long <- with(dat_allo,
             allocate(
               F = F,
               Q = Q,
               alpha = alpha,
               w = w,
               K = K,
               eps_K = .0001,
               eps_lam = .0001,
               Trace = TRUE
             ))
lam_iters <- long$lambdas

dat_allo <- dat_allo %>% mutate(Lambda_adj = pmap(., margexb_fun)) %>%
  mutate(Lambda = pmap(.[names(.) != "q_alpha"], margexb_fun))

K_funfac <- function(df = dat_allo, forecasts = NULL) {
  if (!is.null(forecasts)) df <- df[forecasts,]
  Vectorize(with(df, function(lambda) {
    sum(w * pmax(0, map2_dbl(Q, pmax(0, alpha-w*lambda), exec)))
  }))
}
K_fun <- K_funfac()

ymax = .25
n = nrow(dat_allo)
with(dat_allo,
{
  p1 <- ggplot() + xlim(0, 1) + ylim(0, ymax) +
    map(1:n, ~ geom_function(aes(color = name[[.]]), fun = Lambda_adj[[.]])) +
    geom_hline(yintercept = 0) +
    map(lam_iters, ~ geom_hline(yintercept = ., alpha = .3)) +
    map(1:n, ~ geom_vline(aes(xintercept = (allo / q_alpha)[.], color = name[[.]]), alpha=.5)) +
    labs(color = "Distribution",
         y = expression(lambda(x)),
         x = expression(x[i] / q[list(F[i],alpha[i])])) +
    theme_classic() +
    theme(legend.position = "left")

  p2 <- ggplot(data = tibble(y = seq(0, 1, length.out = 1000)), aes(x = x, y = y)) +
    ylim(0, ymax) + xlim(0, K_fun(0)) +
    geom_line(aes(x = K_funfac()(y))) +
    geom_line(aes(x =  K_funfac(forecasts = c(2, 4))(y)), color = "blue") +
    # #geom_line(aes(x =  K_funfac(forecasts = c(1,2,4))(y))) +
    #map(1:n, ~ geom_line(aes(x = K_funfac(forecasts = c(.))(y), color = name[[.]]))) +
    geom_vline(xintercept = K) +
    labs(x = expression(w ^ T * Q(alpha - w * lambda))) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )+
    theme(legend.position = "none")

p1+p2
})

#tibble(y = seq(0, .2, length.out = 10), x= K_funfac(forecasts = c(1))(y))
