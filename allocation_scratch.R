library(tidyverse)
mus <- 10 + c(20,2,3,4,5)
sigs <- c(2,4,3,1,4)


lmc <- function(Ovg, Und, w, K) {
  Vectorize(function(lambda) {
    sum(qnorm((Und - w*lambda)/(Und + Ovg), mus[1:5], sigs[1:5])) - K}
  )
}

out <- function(Ovg, Und, w, lambda) {
  qnorm((Und - w*lambda)/(Und + Ovg), mus[1:5], sigs[1:5])
}

lmc1 <- lmc(1,1,1,100)
lmc1(2)

uniroot(lmc(10,5,1,85), c(-10,5), extendInt = "yes")

ggplot() + xlim(-10,5) + geom_function(fun = lmc(10,5,1,0), n = 1000)

out(0,1,1,.9)

ggplot() + xlim(0,1) + geom_function(fun = function(x) qnorm(x,1,19), n=1000) +
  geom_function(fun = function(x) qnorm(x,0,10) + qnorm(x,0,3) + qnorm(x,1,6), n=1000, color = "red")

ggplot() + xlim(-10,10) + geom_function(fun = function(x) dnorm(x,2,sqrt(8)), n=1000) +
  geom_function(fun = function(x) dnorm(x,1,2) + dnorm(x,1,2), n=1000, color = "red")

norm_loss <- function(alpha, mean, sd) {
  f <- function(x){
    alpha^(-1)*((pnorm(x, mean, sd)-alpha)*(x-mean) + sd^2*dnorm(x,mean,sd))
  }
  return(f)
}
sig = .2
mag = 10
alphas =   c(.1,.4,.5,.9)
w <- 6

plot_loss <- function(
    alphas,
    lossfun, 
    parms, 
    pfun, 
    qfun,
    mag, 
    xl, 
    xr,
    yb,
    yexpand = c(1,1)) {
  p <- ggplot() + xlim(xl,xr)+ map(
    alphas,
    function(alpha) {
      geom_function(aes(color = as_factor(alpha)), 
                    fun = rlang::exec(lossfun, alpha, !!!parms))
    }
  ) + map(
    alphas,
    function(alpha) {
      geom_segment(aes(color = as_factor(alpha), 
                       x = rlang::exec(qfun, alpha, !!!parms), 
                       y = yb,
                       xend = rlang::exec(qfun, alpha, !!!parms),
                       yend = mag*alpha),
                   alpha = .5)
    }
  ) + map(
    alphas,
    function(alpha) {
      geom_segment(aes(color = as_factor(alpha),
                       x = rlang::exec(qfun, alpha, !!!parms),
                       y = mag*alpha,
                       xend = xl, 
                       yend = mag*alpha),
                   alpha = .5)
    }
  ) + 
    geom_function(fun = function(q) mag*rlang::exec(pfun, q, !!!parms)) +
    scale_y_continuous(expand = yexpand,
                       labels = function(y) y/mag, 
                       breaks = mag*alphas, 
                       limits = c(yb, mag*1.1)) +
    scale_x_continuous(expand = c(0,.1), limits = c(xl,xr)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank())
  return(p)
}



plot_loss(
    alphas = alphas,
    lossfun = norm_loss, 
    parms = list(mean = 0, sd = 2), 
    pfun = pnorm, 
    qfun = qnorm,
    mag = mag, 
    xl = -w, 
    xr = w)


unif_loss <- function(alpha, a, b) {
  f <- function(x){
    (alpha*2*(b-a))^(-1)*(x^2 - 2*(a+(b - a)*alpha)*x+(1 - alpha)*a^2 +alpha*b^2)
  }
  return(f)
}

plot_loss(
  alphas = alphas,
  lossfun = unif_loss, 
  lossparms = list(a = 0, b=1), 
  pfun = punif, 
  qfun = qunif,
  distparms = list(min = 0, max = 1),
  mag = 1.2, 
  xl = -.3, 
  xr = 1.3)

exp_loss <- function(alpha, rate) {
  f <- function(x){
    alpha^(-1)*(exp(-x*rate)*rate^(-1) - (alpha-1) * x - rate^(-1))
  }
  return(f)
}
parms = list(rate = 1)

plot_loss(
  alphas = alphas,
  lossfun = exp_loss, 
  parms = parms, 
  pfun = pexp, 
  qfun = qexp,
  mag = 3, 
  xl = -.5, 
  xr = 5,
  yexpand = c(0,.1),
  yb = -1
  )

unif_loss2 <- function(alpha1,alpha2,a1,a2,b1,b2) {
  f <- function(x1,x2) {
    unif_loss(alpha1, a1, b1)(x1) + unif_loss(alpha2, a2, b2)(x2)
  }
  return(f)
}
x1 <- seq(0, 1, length.out = 100)
x2 <- seq(0, 1, length.out = 100)
alpha1 = .8
alpha2 = .2
a1 = 0
a2 = 0
b1 = 1
b2 = 1
dat <- expand_grid(x1, x2) %>% mutate(z = unif_loss2(alpha1,alpha2,a1,a2,b1,b2)(x1, x2))
ggplot(dat, aes(x1,x2)) + geom_contour(aes(z=z))

norm_loss <- function(alpha, mu, sigma) {
  f <- function(x){
    alpha^(-1)*((pnorm(x, mu, sigma)-alpha)*(x-mu) + sigma^2*dnorm(x,mu,sigma))
  }
  return(f)
}
norm_loss2 <- function(alpha1, mu1, sigma1, alpha2, mu2, sigma2) {
  f <- function(x1,x2) {
    norm_loss(alpha1, mu1, sigma1)(x1) + norm_loss(alpha2, mu2, sigma2)(x2)
  }
}
x1 <- seq(-1, 1, length.out = 200)
x2 <- seq(-1, 1, length.out = 200)
alpha1 = .8
alpha2 = .6
mu1 = .5
mu2 = 0
sigma1 = .5
sigma2 = .1
dat <- expand_grid(x1, x2) %>% 
  mutate(z = unif_loss2(alpha1,alpha2,mu1,mu2,sigma1,sigma2)(x1, x2))
ggplot(dat, aes(x1,x2)) + geom_contour(aes(z=z), bins = 50)
  
# Lau & Lau example in table 1

mu = c(100, 500, 300)
E = c(1,1,2)
U = c(4,1,2)
w = c(1,4,1)


f1 <- function(lambda) pmax(-w*mu*log((1-(U-w*lambda)/(E+U))),0)
f1(2)

f = function(lambda, C) {
  sum(pmax(-w*mu*log((1-(U-w*lambda)/(E+U))),0)) - C
}
lam_star <- uniroot(f, c(0,3))$root

(U-w*lam_star)/(E+U)

Cs = c(1755.2, 1300, 1000, 500, 100, 25)
ggplot() + 
  map(Cs,
      function(C) geom_function(
        fun = Vectorize(function(lam) f(lam, C)), 
        aes(color = factor(C, levels = Cs)))
  )+
  geom_hline(aes(yintercept = 0), alpha = .5)+
  scale_x_continuous(limits = c(-.02,4), 
                     breaks = scales::pretty_breaks(n = 20))+
  scale_color_discrete(breaks = rev(Cs), name = "Constraint")+
  labs(y = "Deficit")+
  theme_bw()

# Lau & Lau example in table 2

a=c(50,300,200)
b=c(150,700,400)

f2_vec <- function(lambda) {a + (b-a)*((U-w*lambda)/(E+U))}
f2_sum <- function(lambda, C) {
  sum(w*f2_vec(lambda)) - C
}
Cs = c(2430, 1600, 400, 200)
ggplot() + 
  map(Cs,
      function(C) geom_function(
        fun = Vectorize(function(lam) f2_sum(lam, C)), 
        aes(color = factor(C, levels = Cs)))
  )+
  geom_hline(aes(yintercept = 0), alpha = .5)+
  scale_x_continuous(limits = c(-.02,2.5), 
                     breaks = scales::pretty_breaks(n = 20))+
  scale_color_discrete(breaks = rev(Cs), name = "Constraint")+
  labs(y = "Deficit")+
  theme_bw()

(2430-400)/(sum((b-a)*w^2/(E+U)))
(2430-1600)/(sum((b-a)*w^2/(E+U)))
alp=U/(E+U)
Z = function(x) ((E+U)/(2*(b-a)))*(x^2-2*(a+(b-a)*alp)*x+(1-alp)*a^2+alp*b^2)
Z(f2_vec(0.2538226)) %>% sum()
Z(c(125,296.9,287.5)) %>% sum()
