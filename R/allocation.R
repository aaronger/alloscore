#' @importFrom purrr map map2 pmap map_dbl map2_dbl map_int partial map_lgl
#' @importFrom rlang exec is_missing is_list caller_env
#' @importFrom dplyr mutate arrange
#' @importFrom tibble tibble
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
#' @param df data frame containing columns for other list arguments which are used only when
#'  those arguments (e.g. `F`) are not passed directly
#' @param F list of cdf functions for forecast distributions
#' @param Q list of quantile functions for forecast distributions
#' @param w numeric vector with cost per unit resource allocated to each
#'  coordinate
#' @param K vector of constraints on total provision
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param dg numeric constant(s) or function(s) to calculate the derivative of the
#'  gpl function `g` for each coordinate
#' @param eps_K
#' @param eps_lam
#' @param point_mass_window
#'
#' @return data frame with (list-)columns
#' \describe{
#'   \item{K}{constraints}
#'   \item{x}{allocations for constraints}
#'   \item{xs}{array with allocations at each iteration as rows}
#'   \item{lam_seq}{sequence of Lagrange multipliers tested at each iteration}
#' }

#' @export
#'
#' @examples
allocate <- function(df = NULL, F, Q, w = 1, K,
                     kappa = 1, alpha = 1,
                     dg = 1,
                     eps_K = .01,
                     eps_lam = 1e-5,
                     point_mass_window = .001,
                     Trace = FALSE,
                     verbose = FALSE) {
  if (!is.null(df)) get_args_from_df(df)
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

  # initialize allocation at q_i if alpha_i < 1 or something big enough to
  # violate constraint if not
  qs <- map2_dbl(Q, alpha, exec)
  qs[qs==Inf] <- (w^(-1)*rep(1,N)*2*max(K))[qs==Inf]

  # set up a data frame to track iterations for all K
  Kdf <- tibble(
    K = K,
    xs = list(qs), # for iterations
    x = xs, # for final allocation
    qs_OK = map(x, ~sum(w*.)) < K,
    converged = FALSE
  ) %>% arrange(K)

  out <- Kdf %>% dplyr::filter(qs_OK) %>% mutate(converged = TRUE)
  Kdf <- Kdf %>% dplyr::filter(!qs_OK)
  # return if quantiles satisfy all constraint levels
  if (nrow(Kdf) == 0) {
    return(out)
  }
  # if not, get lmabda_i'a
  Lambda <- pmap(largs[names(largs) != "Q"], margexb_fun)
  # store ther values at boundary
  Lambda0 <- map2_dbl(Lambda, 0, exec)
  #create columns for binary search interval endpoints at lambda = MEB = 0 (corresponding to
  # allocation at quantile) and the maximum over target set of MEB's evaluated at x = 0
  # (corresponding to the maximum MEB if the allocations were all reduced as much as possible - to 0)
  # as well as columns for the current lamda and previous lambda iterates
  Kdf <- Kdf %>% mutate(
    lamL = 0,
    lamL_seq = list(0),
    lamU = max(Lambda0)*(1+2*eps_lam),
    lamU_seq = list(lamU),
    lam = 0,
    lam_seq = list(0),
    lam_prev = 0
  )
  # store

  # create counter for labeling iterates
  tau <- 1

  # return x = 0 if there is no marginal benefit to allocating more than 0 in any component
  if (Kdf$lamU[1] <= 0) {
    Kdf$x <- list(rep(0, N))
    out <- rbind(out, Kdf) %>% arrange(K) # would out ever be non-empty here?
    return(out)
  }

  while (nrow(Kdf) > 0) {
    Kdf <- Kdf %>% mutate(
      # calculate next lambdas
      lam = (lamL + lamU)/2,
      # store them
      lam_seq = map2(lam_seq, lam, function(lam_seq, lam) c(lam_seq, lam))
    )
    lams <- unique(Kdf$lam)
    for (lam_tau in lams) {
      Kdf_lam <- Kdf %>% dplyr::filter(lam == lam_tau)
      x_tau <- Kdf_lam$x[[1]] # this is the whole point of using Kdf
      stopifnot(length(lam_prev <- unique(Kdf_lam$lam_prev)) == 1) # sanity check
      for (i in 1:N) {
        if (lam_tau > Lambda0[i]) {
          # then lam gives a negative quantile on i so we allocate nothing at this iteration
          x_tau[i] <- 0
        } else {
          # we should expect a critical pt greater than x_tau[i] = 0 and write the ith
          # component of the gradient equation
          lam_grad = function(xi) {
            Lambda[[i]](xi) - lam_tau
          }
          # obtain a search interval in which to look for a root of lam_grad
          if (lam_tau < lam_prev) {
            # if lam has decreased we need to look closer to the quantile
            I <- c(x_tau[i] * (1 - point_mass_window), qs[i])
          } else {
            # and if not we need to look further toward 0
            I <- c(0, x_tau[i] * (1 + point_mass_window))
          }
          if (lam_grad(I[1]) * lam_grad(I[2]) < 0) {
            # then adjust x_tau[i] to that root.
            x_tau[i] <- uniroot(f = lam_grad, interval = I)$root
            # tryCatch(
            #   x_tau[i] <- uniroot(f = lam_grad, interval = I)$root,
            #   error = function(e) {
            #     message("error at tau = ", tau, "  i = ", i)
            #     browser()
            #     stop(e)
            #   }
            # )
          }
          # Otherwise leave x_tau[i] unchanged which will cause lam to increase on next step
        }
      }
      Kdf_lam <- Kdf_lam %>% mutate(
        x = list(x_tau),
        xs = map(xs, ~rbind(., x_tau)),
        lam_prev = lam_tau,
        # check whether we have under- or over-shot K and narrow search interval accordingly
        lamL = ifelse(sum(w*x_tau) > K, lam_tau, lamL),
        lamL_seq = map2(lamL_seq, lamL, function(lamL_seq, lamL) c(lamL_seq, lamL)),
        lamU = ifelse(sum(w*x_tau) <= K, lam_tau, lamU),
        lamU_seq = map2(lamU_seq, lamU, function(lamU_seq, lamU) c(lamU_seq, lamU)),
        # check if converged
        converged = ((lamU - lamL)/(lamU) < eps_lam) | (lamU <= eps_lam)
      )
      Kdf[Kdf$lam == lam_tau,] <- Kdf_lam
    }
    out <- rbind(out, Kdf %>% dplyr::filter(converged))
    Kdf <- Kdf %>% dplyr::filter(!converged)
    tau <- tau + 1
  }

    # old code for dealing with some convergence issue
    # if (tau > 25) {
    #   if (xis %>% map_dbl(~sum(abs(diff(tail(.,3))))) %>% max() < .001) break
    # }
  out <- out %>% mutate(x = pmap(list(x, K, xs),
    function(x, K, xs) {
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
          browser()
          x_L <- x_U <- x
          iter_num <- dim(xs)[1]
          for (j in 1:iter_num) {
            if (j == iter_num + 1) {
              stop("Something went very wrong")
            }
            x_L <- pmin(x_L, xs[[iter_num - j]])
            x_U <- pmax(x_U, xs[[iter_num - j]])
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
      return(x)
    }))
  return(out)
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
                             g = function(u) u,
                             components = FALSE) {
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
  component_scores <- gpl(oracle_allos, y)
  score <- sum(component_scores)
  if (components) {
    return(list(score = score, components = component_scores))
  }
  return(score)
}

#' Obtain the allocation score for a given forecast distribution F for the
#' observed data value y in a constrained allocation problem.
#'
#' @param df data frame containing columns for other list arguments which are used only when
#'  those arguments (e.g. `F`) are not passed directly
#' @param y numeric observed data value
#' @param g list of functions that calculate the gpl function for each coordinate;
#'  default value of NULL causes pinball loss to be used
#' @param against_oracle logical; if `TRUE`, scores are normalized relative to
#'  an oracle forecaster
#' @param components logical; if TRUE, the components of the score (relative to the components
#' of the oracle score if `against_oracle` is `TRUE`) are also returned
#' @inheritParams allocate
#' @export
alloscore <- function(df = NULL, y, F, Q, w = 1, K,
                      kappa = 1, alpha = 1,
                      dg = 1,
                      eps_K = .01,
                      eps_lam = 1e-5,
                      g = NULL,
                      against_oracle = TRUE,
                      verbose = FALSE,
                      components = FALSE) {
  if (!is.null(df)) get_args_from_df(df)
  if (dg == 1) {
    if (!is.null(g)) {
      stop("derivatives of non-identity gpl functions must be specified")
    }
    g <- function(u) u
  }
  allos <- allocate(
    F = F,
    Q = Q,
    w = w,
    K = K,
    kappa = kappa,
    alpha = alpha,
    dg = dg,
    eps_K = eps_K,
    eps_la = eps_lam,
    verbose = verbose)
  gpl <- gpl_loss_fun(g, kappa, alpha)
  component_scores <- gpl(allos, y)
  score <- sum(component_scores)
  if (against_oracle) {
    oscore <- oracle_alloscore(y, w, K,
                               kappa, alpha,
                               dg, eps_K, eps_lam, g,
                               components = components)
    score <- score - oscore[[1]]
    if (components) {
      component_scores <- component_scores - oscore[[2]]
    }
  }
  if (components) {
    return(list(score = score, components = component_scores))
  }
  return(score)
}
