library(tidyverse)
library(patchwork)
library(devtools)
load_all()

(dat <- tibble(
  alpha = c(.2,.6,.3,.9,.5),
  w = c(3,2,3,5,2),
  dists_and_params = list(
    list(dist = "norm", mean = 3, sd = 1),
    list(dist = "norm", mean = 2, sd = 4),
    list(dist = "exp", rate = .5),
    list(dist = "gamma", shape = 1, rate = .3),
    list(dist = "unif", min = 5, max = 6)
    #list(dist = "distfromq", ps = c(.1,.3,.6,.9,.95), qs = c(2,3,3,4,5))
  ),
  name = map(dists_and_params,
             ~toString(paste0(.[["dist"]], ", ",
                              toString(paste(names(.)[-1], .[-1], sep = "=")))))) %>%
  unnest_wider(col = dists_and_params) %>%
  add_pdqr_funs(types = c("p", "q")) %>% mutate(
    q_scale = map2_dbl(Q, alpha, exec),
    # add marginal expected benefit functions scaled to fit on single [0,1] interval,
    Lambda_scaled = pmap(., margexb_fun),
    # and unscaled as well, which are not yet being used
    Lambda = pmap(.[names(.) != "q_scale"], margexb_fun),
    Eloss = pmap(.[c("alpha", "F")], gpl_loss_exp_fun)
  )
)

make_iter_pic(2, dat)

make_iter_pic <- function(K, dat) {
  (dat_allo <- dat %>%
     mutate(
       allo = allocate(
         F = F,
         Q = Q,
         alpha = alpha,
         w = w,
         K = K,
         eps_K = .01,
         eps_lam = .0001
       ),
       allo_q_scaled = allo/q_scale
     ) %>%
     relocate(allo, allo_q_scaled, q_scale)
  )

  # get all iterations of lambda from allocation algorithm
  long <- with(dat_allo,
               allocate(
                 F = F,
                 Q = Q,
                 alpha = alpha,
                 w = w,
                 K = K,
                 eps_K = .01,
                 eps_lam = .0001,
                 Trace = TRUE
               ))
  lam_iters <- long$lambdas


  # make function of lambda that gives total resource use at "lambda-quantiles"
  # which are regular quantiles when lambda = 0.
  # When `!is.null(forecasts))`, the total use is for the subset of rows given by
  # vector `forecasts`
  K_fun_factory <- function(df = dat_allo, forecasts = NULL) {
    if (!is.null(forecasts)) df <- df[forecasts,]
    Vectorize(with(df, function(lambda) {
      alpha_lam <- alpha-w*lambda
      positives <- alpha_lam>0
      qs_plus <- pmax(0, map2_dbl(Q, alpha_lam*positives, exec))
      qs <- qs_plus*positives
      sum(w * qs)
    }))
  }
  K_fun <- K_fun_factory()

  ymax = max(map_dbl(dat_allo$Lambda, ~.(0)))*(1.1)
  xmax = max(dat_allo$allo)*1.1
  n = nrow(dat_allo)
  with(dat_allo,
       {
         p1 <- ggplot() + xlim(c(0,xmax)) + ylim(0, ymax) +
           geom_hline(yintercept = 0) +
           map(lam_iters, ~ geom_hline(yintercept = ., alpha = .3)) +
           geom_label(data = enframe(lam_iters[2:5]), aes(x=0,y=value,label=name)) +
           map(1:n, ~ geom_function(aes(color = name[[.]]), fun = Lambda[[.]])) +
           map(1:n, ~ geom_vline(aes(xintercept = allo[.], color = name[[.]]), alpha=.5)) +
           labs(color = "Distribution",
                y = expression(lambda(x)),
                x = expression(x[i])) +
           theme_classic() +
           theme(legend.position = "left")

         p2 <- ggplot(data = tibble(y = seq(0, 1.01, length.out = 1000), x = K_fun(y)),
                      aes(x = x, y = y)) +
           ylim(0, ymax) + xlim(0, K_fun(0.01)) +
           geom_line() +
           geom_label(data = enframe(lam_iters[2:5]) %>% mutate(x = K_fun(value)), aes(x=x,y=value,label=name)) +
           # geom_line(aes(x =  K_funfac(forecasts = c(2, 4))(y)), color = "blue") +
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
}




# Plot MEBs
# with(dat,
#      {p <- ggplot() + xlim(-1, 10) +
#          map(1:nrow(dat), ~ geom_function(aes(color = name[[.]]), fun = Lambda[[.]]))
#      p}
# )

plotdat <- function(dat) {
  xl = 0
  xr = 10
  plt_dat <- dat
  p <- ggplot(data = data.frame(
    x = seq(from = xl, to = xr, length.out = 1000)
  ))
  p <- p + map(1:nrow(plt_dat), function(i) geom_line(
    aes(x = x, y = plt_dat$F[[i]](x), color = plt_dat$name[[i]]), linetype = 5))
  p_lam <- p + map(1:nrow(plt_dat), function(i) geom_line(
    aes(x = x, y = plt_dat$Lambda[[i]](x), color = plt_dat$name[[i]])))
  #p_lam <- p_lam+ylim(c(-1,1))
  p_loss <- p + map(1:nrow(plt_dat), function(i) geom_line(
    aes(x = x, y = plt_dat$Eloss[[i]](x), color = plt_dat$name[[i]])))
  #p_loss <- p_loss+ylim(c(0,2))
  p_lam + p_loss + plot_layout(ncol=1)
}

plotdat(dat)
