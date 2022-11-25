library(tidyverse)

zxh1 <- read_csv("zxh1.csv") 

zxh <- zxh1 %>% mutate(
  F = pmap(list(mu, sigma), function(mu,sigma) {function(x) {pnorm(x, mean = mu, sd = sigma)}}),
  f = pmap(list(mu, sigma), function(mu,sigma) {function(x) {dnorm(x, mean = mu, sd = sigma)}}),
  Q = pmap(list(mu, sigma), function(mu,sigma) {function(a) {qnorm(a, mean = mu, sd = sigma)}}),
  kaparams = pmap(list(c,h,v), 
                  function(c,h,v) stdize_news_params(ax = c, a_minus = v, a_plus = h))
  ) %>% unnest_wider(kaparams) %>% 
  mutate(
    q = map2_dbl(Q, alpha, exec),
    gpl_exp_loss = pmap(list(kappa, alpha, c, mu, f),
                        function(kappa, alpha, c, mu, f) {
                          gpl_loss_exp_fun(
                            kappa = kappa,
                            alpha = alpha,
                            const = c*mu,
                            f = f)
                        }
    )
  )

zxh <- within(zxh, allo <- allocate(
  F = F,
  Q = Q,
  kappa = kappa,
  alpha = alpha,
  w = c,
  dg = 1,
  K = 2500, 
  eps_K = .001,
  eps_lam = .01
))

zxh_long <- with(zxh, allocate(
  F = F,
  Q = Q,
  kappa = kappa,
  alpha = alpha,
  w = c,
  dg = 1,
  K = 2500, 
  eps_K = .001,
  eps_lam = .01,
  Trace = TRUE
))

zxh <- zxh %>% mutate(
  Z_Opt = map2_dbl(gpl_exp_loss, Opt, exec),
  Z_allo = map2_dbl(gpl_exp_loss, allo, exec)
)

zxh %>% summarise(Z_Opt = sum(Z_Opt), Z_allo = sum(Z_allo))

