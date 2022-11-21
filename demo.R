library(tidyverse)

zxh1 <- read_csv("zxh1.csv") %>% mutate(
  alpha = (v-c)/(v+h),
  F = pmap(list(mu, sigma), function(mu,sigma) {function(x) {pnorm(x, mean = mu, sd = sigma)}}),
  f = pmap(list(mu, sigma), function(mu,sigma) {function(x) {dnorm(x, mean = mu, sd = sigma)}}),
  Q = pmap(list(mu, sigma), function(mu,sigma) {function(a) {qnorm(a, mean = mu, sd = sigma)}}),
  q = map2_dbl(Q, alpha, exec),
  gpl_exp_loss = pmap(list(v, h, c, f, mu),
                      function(v, h, c, f, mu) {
                        gpl_loss_exp_fun(
                          Und = v-c,
                          Ovg = h+c,
                          const = c*mu,
                          f = f)
                      }
  )
)

zxh1$allo <- allocate(
  F = zxh1$F,
  Q = zxh1$Q,
  alpha = zxh1$alpha,
  w = zxh1$c,
  dg = zxh1$v - zxh1$c,
  K = 2500, 
  eps = .00001
)

zxh1_long <- allocate(
  F = zxh1$F,
  Q = zxh1$Q,
  alpha = zxh1$alpha,
  w = zxh1$c,
  dg = zxh1$v - zxh1$c,
  K = 2500, 
  eps = .00001,
  Trace = TRUE
)

zxh1 <- zxh1 %>% mutate(
  Z_Opt = map2_dbl(gpl_exp_loss, Opt, exec),
  Z_allo = map2_dbl(gpl_exp_loss, allo, exec)
)

zxh1 %>% summarise(Z_Opt = sum(Z_Opt), Z_allo = sum(Z_allo))
