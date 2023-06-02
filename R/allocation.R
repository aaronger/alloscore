#' @importFrom purrr map map2 pmap map_dbl map2_dbl map_int partial map_lgl
#' @importFrom rlang exec is_missing is_list caller_env
#' @importFrom dplyr mutate arrange
#' @importFrom tibble tibble
NULL

#' Allocate to minimize expected gpl loss under forecasts F with constraint K
#'
#' @param df data frame containing columns for other list arguments which are used when
#'  those arguments (e.g. `F`) are left empty or give the name of a column of df,
#'  e.g., `F = "cdf"`
#' @param K vector of constraints on total provision (cannot be supplied via df)
#' @target_names vector of names for allocation targets; will default to indices corresponding to
#'  other vector arguments
#' @param F list of cdf functions for forecast distributions
#' @param Q list of quantile functions for forecast distributions
#' @param w numeric vector with cost per unit resource allocated to each
#'  coordinate
#' @param kappa scale factor
#' @param alpha normalized loss when outcome y exceeds forecast x
#' @param g list of strings describing increment functions which will be used
#'  to form gpl functions for each coordinate; should be written using `x`, e.g. "log(x)";
#'  will be differentiated using `stats::D` so don't try anything fancy; defaults to NULL
#'  deferring to `dg`
#' @param dg numeric constant(s) or function(s) to calculate the derivative of the
#'  increment function `g` for each coordinate; defaults to 1 corresponding to pinball loss
#' @param eps_lam
#' @param eps_K
#' @param point_mass_window
#'
#' @return data frame with (list-)columns
#' \describe{
#'   \item{K}{constraints}
#'   \item{x}{allocations for constraints}
#'   \item{xs}{data frames with allocations at each iteration as columns}
#'   \item{lam_seq}{sequences of Lagrange multipliers tested at each iteration}
#'   \item{qs_OK}{}
#'   \item{post-processed}{whether post-processing was used at this K to deal with
#'    plateaus in the objective function; TRUE in particular whenever oracle allocation is
#'    performed for an active K.}
#' }

#' @export
#'
#' @examples
allocate <- function(df = NULL, K,
                     target_names = NA,
                     F, Q, w = 1,
                     kappa = 1, alpha = 1,
                     g = "x", dg = NA,
                     eps_lam = 1e-4,
                     eps_K = .01,
                     point_mass_window = .001) {
  if (is.character(target_names) && (length(target_names) == 1)) {
    target_col_name <- target_names
  } else {
    target_col_name <- NULL
  }
  if (!is.null(df)) {
    get_args_from_df(df)
    N <- nrow(df)
  } else {
    N <- length(F)
  }
  gpl <- new_gpl_df(N = N, target_names = target_names, g = g, kappa = kappa, alpha = alpha)
  if (!is.null(target_col_name)) {
    gpl <- gpl %>% rename(!!target_col_name := target_names)
  } else {
    target_col_name <- "target_names"
  }

  # initialize allocation at q_i if alpha_i < 1 or something big enough to
  # violate constraint if not
  qs <- map2_dbl(Q, alpha, exec)
  qs[qs==Inf] <- (w^(-1)*rep(1,N)*2*max(K))[qs==Inf]
  # get marginal expected benefits
  Lambda <- meb_gpl_df(df = gpl, F = F, w = w)
  # store ther values at boundary
  Lambda0 <- map2_dbl(Lambda, 0, exec)
  # set up a data frame to track iterations for all K
  # with columns for binary search interval endpoints at lambda = MEB = 0 (corresponding to
  # allocation at quantile) and the maximum over target set of MEB's evaluated at x = 0
  # (corresponding to the maximum MEB if the allocations were all reduced as much as possible - to 0)
  # as well as columns for the current lambda and previous lambda iterates
  Kdf <- tibble(
    K = K,
    xs = list(gpl %>% select(target_col_name) %>% mutate(`qs` = qs)), # for iterations
    x = list(tibble::deframe(xs[[1]])), # for final allocation
    qs_OK = map(x, ~sum(w*.)) < K,
    converged = FALSE,
    lamL = 0,
    lamL_seq = list(0),
    lamU = max(Lambda0)*(1+2*eps_lam),
    lamU_seq = list(lamU),
    lam = 0,
    lam_seq = list(0),
    lam_prev = 0
  ) %>% arrange(K)
  # remove rows where quantiles satisfy constraint
  out <- Kdf %>% dplyr::filter(qs_OK) %>% mutate(converged = TRUE)
  Kdf <- Kdf %>% dplyr::filter(!qs_OK)
  # return if there are no other rows
  if (nrow(Kdf) == 0) {
    return(out)
  }
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
        tryCatch({
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
            }
            # Otherwise leave x_tau[i] unchanged which will cause lam to increase on next step
          }
        }, error = function(e) {
          message(paste(
            "Error at index =", i,
            ", target =", target_names[i],
            ", tau =", tau,
            ": ", e$message))
          return(NULL)
        })
      }
      Kdf_lam <- Kdf_lam %>% mutate(
        x = list(x_tau),
        xs = map(xs, ~mutate(., !!as.character(tau) := x_tau)),
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
    out <-tryCatch({rbind(out, Kdf %>% dplyr::filter(converged))},
                   error = function(e) {
                     message(paste("Error at iteration:", tau, "Ks done =", paste(out$K, collapse = ", ")))
                     #browser()
                     })
    Kdf <- Kdf %>% dplyr::filter(!converged)
    tau <- tau + 1
  }
  # post-processing for plateaus, oracles in particular
  out <- out %>% mutate(
    post_processed = pmap_lgl(list(x, K, qs_OK),
                          function(x, K, qs_OK) {
                            (abs((sum(w * x) - K) / K) > eps_K) && !qs_OK
                          }))
  out <- out %>% mutate(
    x = pmap(list(post_processed, x, K, xs, lamU, qs_OK),
             function(post_processed, x, K, xs, lamU, qs_OK) {
               if (post_processed) {
                 if (lamU <= eps_lam) {
                   if (sum(w * x) > K) {
                     stop("strange situation")
                   }
                   Delta <- function(t) {
                     sum(w * t * x) - K
                   }
                   t_star <-
                     uniroot(f = Delta,
                             interval = c(1, 2),
                             extendInt = "upX")$root
                   x <- t_star * x
                 } else {
                   x_L <- x_U <- x
                   iter_num <- ncol(xs) - 2
                   for (j in 1:iter_num) {
                     if (j == iter_num + 1) {
                       stop("Something went very wrong")
                     }
                     x_L <- pmin(x_L, xs[[iter_num - j + 2]])
                     x_U <- pmax(x_U, xs[[iter_num - j + 2]])
                     Delta <- function(t) {
                       sum(w * ((1 - t) * x_L + t * x_U)) - K
                     }
                     if ((Delta(0) < 0) & (Delta(1) > 0))
                       break
                   }
                   t_star <- uniroot(f = Delta, interval = c(0, 1))$root
                   x <- (1 - t_star) * x_L + t_star * x_U
                 }
               }
      return(x)
    }))
  out <- out %>% select(-c(converged, lam_prev)) %>% arrange(K) %>%
    mutate(xdf = map(x, ~tibble::enframe(., name = target_col_name, value = "x")))
  return(structure(
    out,
    class = c("allocated", class(out)),
    gpl_df = gpl,
    w = w,
    target_col_name = target_col_name))
}

#' Score the allocations from set of forecasts against realized outcomes `y`.
#'
#' @param df data frame containing columns for other list arguments which are used only when
#'  those arguments (e.g. `F`) are not passed directly
#'
#' @return a date frame of the form returned by `allocate` with additional columns for
#' \itemize{
#'   \item `components_raw`: the raw gpl losses in each location
#'   \item `score_raw`: the sum of the raw components
#'   \item `components_oracle`: the gpl losses of an oracle in each location
#'   \item `score_oracle`: the sum of the oracle's components
#'   \item `components`: forecaster's difference from the oracle
#'   \item `score`: the forecaster's score minus the oracle's
#' }
#' @export
#'
#' @examples
alloscore <- function(df = NULL, ...) {
  UseMethod("alloscore")
}

#' Default method for alloscore, used when parameters passed individually; works by
#' allocating and then forwarding to allocated method.
#' @param y numeric observed data value
#' @param against_oracle logical; if `TRUE`, components and scores relative to oracle are
#'  included
#' @inheritParams allocate
#' @rdname alloscore
#' @export
alloscore.default <- function(df = NULL, K, target_names = NA,
                              y, F, Q, w = 1,
                              kappa = 1,
                              alpha = 1,
                              g = "x",
                              dg = NA,
                              eps_K = .01,
                              eps_lam = 1e-5,
                              against_oracle = TRUE) {

  # allocate will handle validation and attribute assignments
  if (!is.null(df)) {
    get_args_from_df(df)
  }
  allocate(
    target_names = target_names,
    F = F,
    Q = Q,
    w = w,
    K = K,
    kappa = kappa,
    alpha = alpha,
    g = g,
    dg = dg,
    eps_K = eps_K,
    eps_lam = eps_lam
  ) %>%
    alloscore(y = y, against_oracle = against_oracle) # uses allocated method
}

#' Allocation scoring method for forecasts with computed allocations and gpl loss functions.
#'
#' @param df an allocated data frame
#' @param y numeric observed data value
#' @param against_oracle logical; if `TRUE`, components and scores relative to oracle are
#'  included
#'
#' @rdname alloscore
#' @export
#'
#' @examples
alloscore.allocated <- function(df, y, against_oracle = TRUE) {
  stopifnot(length(y) == length(df$x[[1]]))
  gpl_list <- attr(df, "gpl_df")$gpl_loss_fun # would a join inside the map below be safer?
  scored_df <- df %>% mutate(
    xdf = map(xdf, function(xdf) {
      xdf %>% mutate(
        y = y,
        gpl = gpl_list,
        components_raw = pmap_dbl(list(gpl, x, y), function(gpl, x, y) {gpl(x,y)})
        ) %>% select(-gpl)
    }),
    score_raw = map_dbl(xdf, ~sum(dplyr::pull(.,components_raw)))
  )

  if (against_oracle) {
    # move this code to oracle_alloscore?
    oracle_scores <- attr(df, "gpl_df") %>%
      select(target_names, g, kappa, alpha) %>%
      mutate(kappa = alpha,
             alpha = 1) %>%
      oracle_allocate(y = y,
                      w = attr(df, "w"),
                      K = df$K) %>%
      alloscore(y = y, against_oracle = FALSE) %>%
      mutate(xdf = map(xdf, ~rename(., oracle = x, components_oracle = components_raw))) %>%
      rename(xdf_oracle = xdf, score_oracle = score_raw)
    scored_df <-
      dplyr::left_join(
        scored_df,
        oracle_scores %>% select(K, xdf_oracle, score_oracle), by = "K") %>%
      mutate(
        xdf = map2(xdf, xdf_oracle, function(xdf, xdf_oracle) {
          bind_cols(xdf, xdf_oracle %>% select(oracle, components_oracle)) %>%
            mutate(components = components_raw - components_oracle) %>%
            relocate(oracle, .after = y)
        }),
        score = score_raw - score_oracle
      )
  }
  return(scored_df)
}

#' Allocate according to an oracle's knowledge of outcome y
#'
#' @param y
#' @inheritParams allocate
#'
#' @return see `allocate`
#' @export
#'
#' @examples
oracle_allocate <- function(df = NULL, y, K, ...) {
  allocate(
    df = df,
    F = map(y, function(y) {function(x) {1*(x>=y)}}),
    Q = map(y, function(y) {function(p) {y}}),
    K = K,
    ...
  )
}

#' Score according to an oracle's knowledge of outcome y
#'
#' @param y
#' @inheritParams alloscore
#'
#' @return see `alloscore`
#'
#' @examples
oracle_alloscore <- function(df = NULL, y, K,
                             kappa = 1, alpha = 1,
                             g = function(u) u, ...) {
  gpl <- gpl_loss_fun(g, kappa, alpha)
  oracle_allocate(
    df = df,
    y = y,
    K = K,
    ...
  ) %>% mutate(
    components = map(x, function(x) gpl(x,y)),
    score = map_dbl(components, sum)
  )
}

# Alternative approach to dealing with plateaus in allocate; not currently in use
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
