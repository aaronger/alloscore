#' @importFrom purrr map map2 pmap
#' @importFrom rlang exec is_missing
NULL

#' Convert newsvendor parameters to kappa, alpha
#'
#' @param ax
#' @param ay
#' @param a_minus
#' @param a_plus
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_news_params <- function(ax, ay, a_minus, a_plus) {
  return(list(
    kappa = a_plus + a_minus,
    alpha = (a_minus - ax)/(a_plus + a_minus)
  ))
}

#' Convert over/underprediction parameters to kappa, alpha
#'
#' @param O
#' @param U
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_ou_params <- function(O, U) {
  return(list(
    kappa = O + U,
    alpha = U/(O+U)
  ))
}

#' Convert meteorologist parameters to kappa, alpha
#'
#' @param C
#' @param L
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_met_params <- function(C, L) {
  return(list(
    kappa = L,
    alpha = 1-C/L
  ))
}

#' Create gpl scoring/loss function
#'
#' @param g gpl function
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param U loss when outcome y exceeds forecast x; equals kappa*alpha
#' @param O cost when forecast x exceeds outcome y; equals kappa*(1-alpha)
#' @return function giving loss
#' @export
gpl_loss_fun <- function(
    g = function(u) u,
    kappa = 1,
    alpha,
    O,
    U,
    const = 0) {
  if (!xor(is_missing(U), is_missing(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is_missing(U)) {
    return(
      function(x, y) {
        scale*(O*pmax(g(x) - g(y), 0) + U*pmax(g(y) - g(x),0)) + const
      }
    )
  }
  if (!is_missing(alpha)) {
    return(
      function(x, y) {
        kappa*((1-alpha)*pmax(g(x) - g(y), 0) + alpha*pmax(g(y) - g(x),0)) + const
      }
    )
  }
}

#' Create gpl expected loss function
#'
#' @param ... passed to gpl_loss_fun
#' @param gpl_loss specified loss function if ... missing
#' @param f probability density
#' @export
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

#' Allocate to minimize expected gpl loss under forecasts F with constraint K
#'
#' @param F
#' @param Q
#' @param kappa
#' @param alpha
#' @param dg
#' @param w
#' @param K
#' @param eps_K
#' @param eps_lam
#' @param Trace
#'
#' @return
#' @export
#'
#' @examples
allocate <- function(F, Q, kappa, alpha, dg, w, K, eps_K, eps_lam, Trace = FALSE) {
  # validation
  largs <- list(F = F, Q = Q, kappa = kappa, alpha = alpha, dg = dg, w = w)
  N <-  max(map_int(largs, length))
  for (i in 1:length(largs)) {
    l <- largs[[i]]
    if (length(l) == 1) {
      if (!is_list(l)) {l <- list(l)}
      largs[[i]] <- rep(l, N)
    }
  }
  if (any(!map_int(largs, length) == N)) {
    stop("parameter lists do not have matching lengths")
  }
  # convert constants in dg to functions
  largs$dg <- map(largs$dg, function(dg) {if (is_function(dg)) dg else function(x) dg})
  # make lambda_i's
  Lambda <- pmap(largs[names(largs) != "Q"], function(F, kappa, alpha, w, dg) {
    function(xi) w^(-1)*kappa*dg(xi)*(alpha - F(xi))
  })
  # get quantiles or finite constraint violators if alpha_i = 1
  Q <- map2(Q, w, function(Q, w) {function(alpha) {min(Q(alpha), w^(-1)*2*K)}})
  x <- qs <- map2_dbl(Q, alpha, exec)
  if (Trace) {
    xs <- list(x)
  }
  # return quantiles if they satisfy constraint
  if (sum(w*x) < K) {
    if (Trace) {
      return(list(x = x, xs = xs, lambdas = 0))
    }
    return(x)
  }
  lamL <- 0
  lamU <- max(map2_dbl(Lambda, 0, exec))
  # return 0 if no marginal benefit for any component
  if (lamU <= 0) {
    if (Trace) {
      return(list(x = rep(0, N), xs = xs, lambdas = NA))
    }
    return(rep(0, N))
  }
  lam <- c(0)
  tau <- 2
  # main loop
  while ((abs((sum(w*x) - K)/K) > eps_K) | (lamU - lamL > eps_lam)) {
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


