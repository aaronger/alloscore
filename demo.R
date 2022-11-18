library(tidyverse)

zxh1 <- read_csv("zxh1.csv")
zxh1 <- zxh1 %>% mutate(
  alpha = (v-c)/(v+h),
  F = pmap(list(mu, sigma), function(mu,sigma) {function(x) {pnorm(x, mean = mu, sd = sigma)}}),
  Q = pmap(list(mu, sigma), function(mu,sigma) {function(a) {qnorm(a, mean = mu, sd = sigma)}}),
  q = map2_dbl(Q, alpha, exec)
  )

zxh1_sol <- allocate(
  F = zxh1$F,
  Q = zxh1$Q,
  alpha = zxh1$alpha,
  w = zxh1$c,
  dg = rep(1,17),
  K = 2500, 
  eps = .00001
)
