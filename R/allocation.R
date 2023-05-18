#' @importFrom purrr map map2 pmap map_dbl map_int partial
#' @importFrom rlang exec is_missing is_list
#' @importFrom dplyr mutate
NULL

#' Basic g-linear loss for under-prediction of an outcome y
#'
#' @param g a non-decreasing function
#'
#' @return function of prediction x and outcome y that only penalizes under-prediction
#' @export
#'
#' @examples
under_loss <- function(g = function(u) {u}) {
  function(x,y) {pmax(g(y) - g(x), 0)}
}

#' Basic expected g-linear loss for under-prediction of a random outcome Y
#'
#' @param dg derivative of g, probably will work even with discontinuities
#' @param F predictive cdf of an outcome Y
#' @return function of prediction x giving expected g-linear loss for under-predicting Y~F
#' @export
#'
#' @examples
expected_under_loss <- function(dg = function(y) {1}, F) {
  Vectorize(function(x) {
    integrate(f = function(y) {(1-F(y))*dg(y)}, lower = x, upper = Inf, rel.tol = .001)$value
    })
}

#' Basic g-linear loss for over-prediction of an outcome y
#'
#' @param g a non-decreasing function
#'
#' @return function of prediction x and outcome y that only penalizes over-prediction
#' @export
#'
#' @examples
over_loss <- function(g = function(u) {u}) {
  function(x,y) {pmax(g(x) - g(y), 0)}
}

#' Basic expected g-linear loss for over-prediction of a random outcome Y
#'
#' @param dg derivative of g, probably will work even with discontinuities
#' @param F predictive cdf of an outcome Y
#' @return function of prediction x giving expected g-linear loss for over-predicting Y~F
#' @export
#'
#' @examples
expected_over_loss <- function(dg = function(y) {1}, F) {
  Vectorize(function(x) {
    integrate(f = function(y) {F(y)*dg(y)}, lower = Inf, upper = x, rel.tol = .001)$value
    })
}

#' Create generalized piecewise linear (gpl) scoring/loss function,
#' which in general need not be piecewise linear
#'
#' @param g gpl function
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param U loss when outcome y exceeds forecast x; equals kappa*alpha
#' @param O cost when forecast x exceeds outcome y; equals kappa*(1-alpha)
#' @param const an added constant, defaults to 0 for which gpl loss function `L` will have `L(x,x) = 0`
#' @return function with arguments `x` and `y` giving loss
#' @export
gpl_loss_fun <- function(g = function(u) {u},
                         kappa = 1, alpha, O, U, const = 0) {
  if (!xor(is_missing(U), is_missing(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is_missing(U)) {
    function(x, y) {O*over_loss(g)(x,y) + U*under_loss(g)(x,y) + const}
  }
  if (!is_missing(alpha)) {
    function(x, y) {kappa*((1-alpha)*over_loss(g)(x,y) + alpha*under_loss(g)(x,y)) + const}
  }
}

#' Create gpl expected loss function
#'
#' @param dg derivative of gpl function
#' @param F predictive CDF
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param U loss when outcome y exceeds forecast x; equals kappa*alpha
#' @param O cost when forecast x exceeds outcome y; equals kappa*(1-alpha)
#' @return function with argument x giving the expected loss with respect to the
#'  distribution F
#' @export
gpl_loss_exp_fun <- function(dg = function(u) {1}, F,
                              kappa = 1, alpha, O, U,
                              const = 0) {
  if (!xor(is_missing(U), is_missing(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is_missing(U)) {
    function(x) {O*expected_over_loss(dg,F)(x) + U*expected_under_loss(dg,F)(x)}
  }
  if (!is_missing(alpha)) {
    function(x) {kappa*((1-alpha)*expected_over_loss(dg,F)(x) + alpha*expected_under_loss(dg,F)(x))}
  }
}

#' Create function to calculate the marginal expected benefit of allocating an
#' additional unit of resources to a given target.
#'
#' @param q_scale scaling factor for resource level x; used for aligning marginal expected
#'  benefit functions of different targets (and parameters) on single x interval.
#' @param ... ignored
#' @inheritParams allocate
#'
#' @return a function with argument x that calculates the marginal expected benefit
#' @export
#' @details
#' Note that this is minus the derivative of the expected score \eqn{\overline{s}_F}
#'
#' @examples
margexb_fun <- function(F, kappa = 1, alpha, w, dg = 1, q_scale = 1,...) {
  if (!is_function(dg)) {
    return(function(x) {w^(-1)*kappa*dg*(alpha - F(x*q_scale))})
  }
  function(x) {w^(-1)*kappa*dg(x*q_scale)*(alpha - F(x*q_scale))}
}

find_linear_intervals <- function(Lambda, q, grid_size = 1000, tol = min(.001, 1/grid_size)) {
  intervals <- list()
  k <- 1
  dtol <- tol*abs(Lambda(0)-Lambda(q))
  x <- seq(from = 0, to = q, length.out = grid_size)
  i <- 2
  while (i < grid_size) {
    if (abs(Lambda(x[i])-Lambda(x[i-1])) < dtol) {
      left <- x[i-1] # start interval
      i <- i+1 # give it at least length one
      while ((abs(Lambda(x[i])-Lambda(left)) < dtol) & (i < grid_size)) {
        i <- i+1 # keep adding if we don't leave dtol nbd
      }
      right <- x[i-1] # end interval at last point in dtol nbd
      lam_val <- mean(c(Lambda(left), Lambda(right)))
      intervals[[k]] <- c(left = left, right = right, lam_val = lam_val)
      k <- k+1
    }
    i <- i+1
  }
  return(intervals)
}

#' Allocate to minimize expected gpl loss under forecasts F with constraint K
#'
#' @param F list of cdf functions for forecast distributions
#' @param Q list of quantile functions for forecast distributions
#' @param w numeric vector with cost per unit resource allocated to each
#'  coordinate
#' @param K constraint on total provision
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param dg numeric constant(s) or function(s) to calculate the derivative of the
#'  gpl function `g` for each coordinate
#' @param eps_K
#' @param eps_lam
#' @param point_mass_window
#' @param Trace logical, record iterations
#' @param Verbose logical, echo K and report number of iterations required to find allocation
#'
#' @return If `Trace` is `FALSE`, a numeric vector of length `N` givig the
#'  optimal allocations. If `Trace` is `TRUE`, a named list with entries:
#'    - `x`: numeric vector of length `N` giving the optimal allocations
#'    - `xs`: list of values of `x` over the course of optimization
#'    - `lambdas`: numeric vector of values of lambda over the course of optimization
#' @export
#'
#' @examples
allocate <- function(F, Q, w, K,
                     kappa = 1, alpha,
                     dg = 1,
                     eps_K, eps_lam,
                     point_mass_window = .001,
                     Trace = FALSE,
                     verbose = FALSE) {
  # validation
  largs <- list(F = F, Q = Q, kappa = kappa, alpha = alpha, dg = dg, w = w)
  N <-  max(map_int(largs, length))
  for (i in 1:length(largs)) {
    l <- largs[[i]]
    # allow for repeated parameters (such as dg = 1) to be given by single term
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
  # initialize allocation at q_i if alpha_i < 1 or something big enough to
  # violate constraint if not
  qs <- map2_dbl(Q, alpha, exec) 
  qs[qs==Inf] <- (w^(-1)*rep(1,N)*2*K)[qs==Inf]
  x <- qs
  xs <- list(x)
  xis <- as.list(x)
  # return this initial allocation at quantiles if it satisfies constraint
  if (sum(w*x) < K) {
    if (Trace) {
      return(list(x = x, xs = xs, lambdas = 0))
    }
    return(x)
  }
  # if not, initialize binary search interval endpoints at lambda = MEB = 0 (corresponding to
  # allocation at quantile) and the maximum over target set of MEB's evaluated at x = 0
  # (corresponding to the maximum MEB if the allocations were all reduced as much as possible - to 0)
  lamL <- 0
  lamU <- max(map2_dbl(Lambda, 0, exec))
  lamU <- lamU*(1+2*eps_lam)
  # return x = 0 if there is no marginal benefit to allocating more than 0 in any component
  if (lamU <= 0) {
    if (Trace) {
      return(list(x = rep(0, N), xs = xs, lambdas = NA))
    }
    return(rep(0, N))
  }
  # initialize sequence of search iterates at the 1st (and now rejected) value 0
  lam <- c(0)
  # create counter for labeling iterates, starting with the soon to be calculated 2nd
  tau <- 1
  # main loop
  while (((lamU - lamL)/(lamU) > eps_lam) & (lamU > eps_lam)) {
    tau <- tau + 1
    lam[tau] <- (lamL + lamU)/2
    for (i in 1:N) {
      if (lam[tau] > Lambda[[i]](0)) {
        # then lam gives a negative quantile on i so we allocate nothing at this iteration
        x[i] <- 0
      } else {
        # we should expect a critical pt greater than x_i = 0 and write the ith
        # component of the gradient equation
        lam_grad = function(xi) {
          Lambda[[i]](xi) - lam[tau]
        }
        # obtain a search interval in which to look for a root of lam_grad
        if (lam[tau] < lam[tau - 1]) {
          # if lam has decreased we need to look closer to the quantile
          I <- c(x[i] * (1-point_mass_window), qs[i])
        } else {
          # and if not we need to look further toward 0
          I <- c(0, x[i] * (1+point_mass_window))
        }
        if (lam_grad(I[1]) * lam_grad(I[2]) < 0) {
        # then adjust x[i] to that root.
          tryCatch(
            x[i] <- uniroot(f = lam_grad, interval = I)$root,
            error = function(e) {
              message("error at tau = ", tau, "  i = ", i)
              browser()
              stop(e)
            }
          )
        }
        # Otherwise leave x[i] unchanged which will cause lam to increase on next step
      }
      xis[[i]][tau] <- x[i]
    }
    if (sum(w*x) < K) {
      lamU <- lam[tau]
    }
    else {
      lamL <- lam[tau]
    }
    xs[[tau]] <- x
    # if (tau > 25) {
    #   if (xis %>% map_dbl(~sum(abs(diff(tail(.,3))))) %>% max() < .001) break
    # }
  }
  if ((abs((sum(w*x) - K)/K) > eps_K)) {
    if (lamU <= eps_lam) {
      if (sum(w*x) > K) {
        stop("strange situation")
      }
      Delta <- function(t) {
        sum(w*t*x) - K
      }
      t_star <- uniroot(f = Delta, interval = c(1,2), extendInt = "upX")$root
      x <- t_star*x
    } else {
      x_L <- x_U <- x
      for (j in 1:tau) {
        if (j == tau + 1) {
          stop("Something went very wrong")
        }
        x_L <- pmin(x_L, xs[[tau - j]])
        x_U <- pmax(x_U, xs[[tau - j]])
        Delta <- function(t) {
          sum(w * ((1 - t) * x_L + t * x_U)) - K
        }
        if ((Delta(0) < 0) & (Delta(1) > 0))
          break
      }
      t_star <- uniroot(f = Delta, interval = c(0, 1))$root
      x <- (1-t_star)*x_L + t_star*x_U
    }
  }
  if (Trace) {
    return(list(x = x, xs = xs, xis = xis, lambdas = lam, meb = Lambda))
  }
  if (verbose) cat(K, tau, " ")
  return(x)
}

oracle_allocate <- function(y, w, K,
                            kappa = 1, alpha,
                            dg = 1,
                            eps_K, eps_lam = .001, Trace = FALSE) {
  allocate(
    F = map(y, function(y) {function(x) {1*(x>=y)}}),
    Q = map(y, function(y) {function(p) {y}}),
    w = w,
    K = K,
    kappa = kappa,
    alpha = alpha,
    dg = dg,
    eps_K = eps_K,
    eps_lam = eps_lam,
    Trace = Trace
  )
}

oracle_alloscore <- function(y, w, K,
                             kappa = 1, alpha,
                             dg = 1,
                             eps_K,
                             eps_lam = .001,
                             g = function(u) u) {
  oracle_allos <- oracle_allocate(
    y = y,
    w = w,
    K = K,
    kappa = kappa,
    alpha = alpha,
    dg = dg,
    eps_K = eps_K,
    eps_lam = eps_lam)
  gpl <- gpl_loss_fun(g, kappa, alpha)
  score <- sum(gpl(oracle_allos, y))
  return(score)
}

#' Obtain the allocation score for a given forecast distribution F for the
#' observed data value y in a constrained allocation problem.
#'
#' @param y numeric observed data value
#' @param g list of functions that calculate the gpl function for each coordinate;
#'  default value of NULL causes pinball loss to be used
#' @param against_oracle logical; if `TRUE`, scores are normalized relative to
#'  an oracle forecaster
#' @inheritParams allocate
alloscore <- function(y, F, Q, w, K,
                      kappa = 1, alpha,
                      dg = 1,
                      eps_K,
                      eps_lam,
                      g = NULL,
                      against_oracle = TRUE,
                      verbose = FALSE) {
  if (dg == 1) {
    if (!is.null(g)) {
      stop("derivatives of non-identity gpl functions must be specified")
    }
    g <- function(u) u
  }
  allos <- allocate(F, Q, w, K,
                    kappa, alpha,
                    dg, eps_K, eps_lam, verbose = verbose)
  gpl <- gpl_loss_fun(g, kappa, alpha)
  score <- sum(gpl(allos, y))
  if (against_oracle) {
    score <- score - oracle_alloscore(y, w, K,
                                      kappa, alpha,
                                      dg, eps_K, eps_lam, g)
  }
  return(score)
}
