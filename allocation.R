library(purrr)
library(rlang)

allocate <- function(F, Q, alpha, dg, w, K, eps, Trace = FALSE) {
# validation
  N <- length(F)
# convert constants in dg to functions
  dg <- map(dg, function(dg) {if (is_function(dg)) dg else function(x) dg})
# make lambda_i's
  Lambda <- pmap(list(F, alpha, w, dg), function(.F, .alpha, .w, .dg) {
    function(xi) .w^(-1)*.dg(xi)*(1-.F(xi)/.alpha)
  })
# get quantiles
  x <- qs <- map2_dbl(Q, alpha, exec)
  if (Trace) {
    xs <- list(x)
  }
  lamL <- 0
  lamU <- max(map2_dbl(Lambda, 0, exec))
  lam <- c(0)
  tau <- 2
# main loop
  while (abs((sum(w*x) - K)/K) > eps) {
    lam[tau] <- (lamL + lamU)/2
    for (i in 1:N) {
      if (lam[tau] <= Lambda[[i]](0)) {
        I <- if (lam[tau] < lam[tau - 1]) c(x[i], qs[i]) else c(0, x[i])
        x[i] <- uniroot(
          f = function(xi) {
            Lambda[[i]](xi) - lam[tau]
          },
          interval = I
        )$root
      }
      else {
        x[i] <- 0
      }
    }
    if (sum(w*x) < K) {
      lamU <- lam[tau]
    }
    else {
      lamL <- lam[tau]
    }
    if (Trace) {
      xs[[tau]] <- x
    }
    tau <- tau + 1
  }
  if (Trace) {
    return(list(x = x, xs = xs, lambdas = lam))
  }
  return(x)
}

#' @param g gpl function
#' @param Und loss (L) when outcome y exceeds forecast x; equals alpha*(Und + Ovg)
#' @param Ovg cost (C) when forecast x exceeds outcome y; equals (1-alpha)*(Und + Ovg)
#' @return function giving loss 
gpl_loss_fun <- function(
    g = function(u) u, 
    Und, 
    Ovg,
    alpha,
    scale = 1,
    const = 0) {
  if (!xor(is_missing(Und), is_missing(alpha))) {
    stop("Either Und or alpha must be specified, but not both")
  }
  if (!is_missing(Und)) {
    if (is_missing(Ovg)) {
      Ovg <- 1 - Und
    }
    return(
      function(x, y) {
        scale*(Ovg*pmax(g(x) - g(y), 0) + Und*pmax(g(y) - g(x),0)) + const
      }
    )
  }
  if (!is_missing(alpha)) {
    return(
      function(x, y) {
        scale*((1-alpha)*pmax(g(x) - g(y), 0) + alpha*pmax(g(y) - g(x),0)) + const
      }
    )
  }
}

gpl_loss_exp_fun <- function(..., gpl_loss = NULL, f) {
  if (is.null(gpl_loss)) {
    gpl_loss <- gpl_loss_fun(...)
  }
  Z <- function(x) {
    map_dbl(x, function(.x)
      integrate(
        f = function(y) gpl_loss(.x, y) * f(y), lower = -Inf, upper = Inf
      )$value)
  }
  return(Z)
}


