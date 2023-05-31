
#' Basic g-linear loss for under-prediction of an outcome y
#'
#' @param g a non-decreasing function
#'
#' @return function of prediction x and outcome y that only penalizes under-prediction
#' @export
#'
#' @examples
under_loss <- function(g = function(u) {u}) {
  if (is.character(g)) {g <- get_function(g)}
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
exp_under_loss <- function(dg = function(y) {1}, F) {
  Vectorize(function(x) {
    integrate(f = function(y) {(1-F(y))*dg(y)}, lower = x, upper = Inf, rel.tol = .001)$value
  })
}

#' Create derivative of expected under-prediction loss
#'
#' @param g
#' @param dg
#'
#' @return
#' @export
#'
#' @examples
dexp_under_loss <- function(g = "x", dg = NULL, F) {
  if (is.null(dg)) {dg <- get_derivative(g)}
  function(x) {dg(x) * (F(x) - 1)}
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
  if (is.character(g)) {g <- get_function(g)}
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
exp_over_loss <- function(dg = function(y) {1}, F) {
  Vectorize(function(x) {
    integrate(f = function(y) {F(y)*dg(y)}, lower = Inf, upper = x, rel.tol = .001)$value
  })
}

#' Create derivative of expected over-prediction loss
#'
#' @param g
#' @param dg
#'
#' @return
#' @export
#'
#' @examples
dexp_over_loss <- function(g = "x", dg = NULL, F) {
  if (is.null(dg)) {dg <- get_derivative(g)}
  function(x) {dg(x) * F(x)}
}

#' Create generalized piecewise linear (gpl) scoring/loss function,
#' which in general need not be piecewise linear
#'
#' @param g a non-decreasing increment function
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param U loss when outcome y exceeds forecast x; equals kappa*alpha
#' @param O cost when forecast x exceeds outcome y; equals kappa*(1-alpha)
#' @param offset an added function of y, defaults to 0 for which gpl loss function `L` will have `L(x,x) = 0`
#' @return function with arguments `x` and `y` giving loss
#' @export
gpl_loss_fun <- function(g = "x", kappa = 1, alpha, O, U = NA, offset = 0) {
  if (is.character(g)) {g <- get_function(g)}
  if (!xor(is.na(U), is.na(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is.na(U)) {
    gpl_base <- function(x, y) {O * over_loss(g)(x,y) + U * under_loss(g)(x,y)}
  }
  if (!is.na(alpha)) {
    gpl_base <- function(x, y) {kappa * ((1-alpha) * over_loss(g)(x,y) +
      alpha * under_loss(g)(x,y))}
  }
  if (is_function(offset)) {
    gpl <- function(x,y) {gpl_base(x,y) + offset(y)}
  } else {
    gpl <- function(x,y) {gpl_base(x,y) + offset}
  }
  return(gpl)
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
exp_gpl_loss_fun <- function(dg = function(u) {1}, F,
                             kappa = 1, alpha, O, U = NA,
                             const = 0) {
  if (!xor(is.na(U), is.na(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is.na(U)) {
    function(x) {O*exp_over_loss(dg,F)(x) + U*exp_under_loss(dg,F)(x)}
  }
  if (!is.na(alpha)) {
    function(x) {kappa*((1-alpha)*exp_over_loss(dg,F)(x) + alpha*exp_under_loss(dg,F)(x))}
  }
}

#' Create derivative of expected gpl loss
#'
#' @param g
#' @param dg
#' @param kappa
#' @param alpha
#' @param O
#' @param U
#'
#' @return function with argument x giving the derivative of expected loss with respect to the
#'  distribution F
#' @export
#'
#' @examples
dexp_gpl_loss <- function(g = "x", dg = NULL, F, kappa = 1, alpha = NA, O, U = NA) {
  if (is_missing(F)) {
    stop("cdf F is missing")
  }
  if (is.null(dg)) {dg <- get_derivative(g)}
  if (!xor(is.na(U), is.na(alpha))) {
    stop("Either U or alpha must be specified, but not both")
  }
  if (!is.na(U)) {
    if (O == 0) { # avoid 0*1/0 problems with g = log(x)
      function(x) {
        U * dexp_under_loss(dg = dg, F = F)(x)
      }
    } else {
      function(x) {
        O * dexp_over_loss(dg = dg, F = F)(x) +
        U * dexp_under_loss(dg = dg, F = F)(x)
      }
    }
  }
  if (!is.na(alpha)) {
    if (alpha == 1) { # avoid 0*1/0 problems with g = log(x)
      function(x) {
        kappa * (alpha * dexp_under_loss(dg = dg, F = F)(x))
      }
    } else {
      function(x) {
        kappa * ((1 - alpha) * dexp_over_loss(dg = dg, F = F)(x) +
                   alpha * dexp_under_loss(dg = dg, F = F)(x))
      }
    }
  }
}

#' Create data frame of gpl functions and associated parameters
#'
#' @param N
#' @param g
#' @param kappa
#' @param alpha
#' @param O
#' @param U
#' @param offset
#'
#' @return
#' @export
#'
#' @examples
new_gpl_df <- function(N = NULL, target_names = NA,
                       g = "x", kappa = 1, alpha = 1, O = NA, U = NA, offset = 0) {
  args <- list(g = g, target_names = target_names,
               kappa = kappa, alpha = alpha, O = O, U = U, offset = offset)
  names <- names(args)
  if (is.null(N)) {
    N <- max(map_int(largs, length))
  }
  for (i in seq_along(args)) {
    if (!length(args[[i]]) %in% c(1, N)) {
      stop(paste(names[i], "must be of length 1 or", N))
    }
  }
  tibble_args <- map(names, function(nm) {
    if (is.list(args[[nm]])) {
      if (length(args[[nm]]) == 1) {
        rep(args[[nm]], N)
      } else {
        args[[nm]]
      }
    } else {
      if (length(args[[nm]]) == N) {
        args[[nm]]
      } else {
        rep(args[[nm]], N)
      }
    }
  })
  names(tibble_args) <- names
  gpl <- as_tibble(tibble_args)
  gpl <- gpl %>% mutate(
    gpl_loss_fun = pmap(list(g, kappa, alpha, O, U, offset), gpl_loss_fun)
  )
  return(structure(gpl, class = c("gpl_df", class(gpl))))
}

#' Create a list of derivatives of expected gpl losses from a gpl_df
#'
#' @param df
#' @param F list of predictive cdfs
#'
#' @return
#' @export
#'
#' @examples
dexp_gpl_df <- function(df, F) {
  stopifnot("gpl_df" %in% class(df))
  if (nrow(df) != length(F)) {
    stop("Need same number of targets and distributions")
  }
  df %>% select(any_of(names(formals(dexp_gpl_loss)))) %>%
    mutate(F = F) %>% pmap(dexp_gpl_loss)
}

#' Create marginal expected benefit functions for vector of gpl losses with respect
#' to predictive distributions F and weights w
#'
#' @param w weights for optimizing allocation
#' @inheritParams dexp_gpl_df
#'
#' @return
#' @export
#'
#' @examples
meb_gpl_df <- function(df, F, w) {
  stopifnot("gpl_df" %in% class(df))
  stopifnot(is.numeric(w))
  if (!(length(w) %in% c(1, nrow(df)))) {
    stop("Need same number of targets and weights")
  }
  df %>% mutate(kappa = -(1/w)*kappa) %>% dexp_gpl_df(F = F)
}





