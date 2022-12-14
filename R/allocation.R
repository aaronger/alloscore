#' @importFrom purrr map map2 pmap partial
#' @importFrom rlang exec is_missing
#' @importFrom dplyr mutate
NULL

#' Create gpl scoring/loss function
#'
#' @param g gpl function
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param U loss when outcome y exceeds forecast x; equals kappa*alpha
#' @param O cost when forecast x exceeds outcome y; equals kappa*(1-alpha)
#' @return function giving loss
#' @export
gpl_loss_fun <- function(g = function(u) {u}, kappa = 1, alpha, O, U, const = 0) {
  if (!xor(is_missing(U), is_missing(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is_missing(U)) {
    gpl_loss <- function(x, y) {
      scale * (O * pmax(g(x) - g(y), 0) + U * pmax(g(y) - g(x), 0)) + const
    }
  }
  if (!is_missing(alpha)) {
    gpl_loss <- function(x, y) {
      kappa * ((1 - alpha) * pmax(g(x) - g(y), 0) + alpha * pmax(g(y) - g(x), 0)) + const
    }
  }
  return(gpl_loss)
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

#' Marginal expected benefit
#'
#' @param kappa
#' @param alpha
#' @param w
#' @param dg
#'
#' @return
#' @export
#'
#' @examples
margexb_fun <- function(F, kappa = 1, alpha, w, dg = 1, q_alpha = 1,...) {
  if (!is_function(dg)) {
    return(function(x) {w^(-1)*kappa*dg*(alpha - F(x*q_alpha))})
    }
  function(x) {w^(-1)*kappa*dg(x*q_alpha)*(alpha - F(x*q_alpha))}
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
#' @param Trace logical, record iterations
#'
#' @return
#' @export
#'
#' @examples
allocate <- function(F, Q, w, K,
                     kappa = 1, alpha,
                     dg = 1,
                     eps_K, eps_lam, Trace = FALSE) {
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
  # get lambda_i's
  Lambda <- pmap(largs[names(largs) != "Q"], margexb_fun)
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

alloscore <- function(y, F, Q, w, K,
                      kappa = 1, alpha,
                      dg = 1,
                      eps_K, eps_lam,
                      g,
                      against_oracle = TRUE) {
  allos <- allocate(F = F, Q = Q, w = w, K = K,
                    kappa = kappa, alpha = alpha,
                    dg = dg,
                    eps_K = eps_K, eps_lam = eps_lam)
  gpl <- gpl_loss_fun(g = g, kappa = kappa, alpha = alpha)
  score <- gpl(allos, y)
}
