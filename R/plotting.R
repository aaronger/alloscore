
plot_loss <- function(
    p = NULL,
    alphas,
    lossfun,
    parms,
    pfun,
    qfun,
    mag,
    xl,
    xr,
    yb,
    yexpand = c(0,.1)) {
  if (is.null(p)) {
    p <- ggplot() + xlim(xl,xr)
  }
  p <- p + map(
    alphas,
    function(alpha) {
      geom_function(aes(color = as_factor(alpha)),
                    fun = rlang::exec(lossfun, alpha, !!!parms), n = 500)
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
    geom_function(fun = function(q) mag*rlang::exec(pfun, q, !!!parms), n=500) +
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


ggplot() + map(datu$F, function(F) geom_function(fun = F, n = 1000)) + xlim(c(-10,15)) +
  map(datu$Eloss, function(Eloss) geom_function(fun = Eloss, n = 1000))

ggplot() + map(dat$f, function(f) geom_function(fun = f, n = 1000)) + xlim(c(-10,15))

xl = 0
xr = 10
plt_dat <- datu
p <- ggplot(data = data.frame(
  x = seq(from = xl, to = xr, length.out = 1000)
))
p <- p + map(1:nrow(plt_dat), function(i) geom_line(
  aes(x = x, y = plt_dat$F[[i]](x), color = plt_dat$name[[i]]), linetype = 5))
p_lam <- p + map(1:nrow(plt_dat), function(i) geom_line(
  aes(x = x, y = plt_dat$Lambda[[i]](x), color = plt_dat$name[[i]])))
p_lam <- p_lam+ylim(c(-1,1))
p_lam
p_loss <- p + map(1:nrow(plt_dat), function(i) geom_line(
  aes(x = x, y = plt_dat$Eloss[[i]](x), color = plt_dat$name[[i]])))
p_loss <- p_loss+ylim(c(0,2))
p_lam + p_loss + plot_layout(ncol=1)


