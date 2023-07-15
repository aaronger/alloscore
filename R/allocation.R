#' @importFrom purrr map map2 pmap map_dbl map2_dbl pmap_dbl map_int partial map_lgl pmap_lgl
#' @importFrom rlang exec is_missing is_list is_function caller_env
#' @importFrom dplyr mutate arrange select filter
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
#' @param eps_lam tolerance for ending iteration on lambda
#' @param eps_K tolerance for post-processing after lambda iteration concludes
#' @param point_mass_window factor to widen search intervals by in root finding so as to
#'  catch point-masses in the F's, intended or otherwise
#'
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
                     F = NULL, Q = NULL, w = 1,
                     kappa = 1, alpha = 1,
                     g = "x", dg = NA,
                     eps_lam = 1e-4,
                     eps_K = .01,
                     point_mass_window = .001) {
  if (is.character(target_names) && (length(target_names) == 1)) {
    # if target_names is a single string meant to identify a column in df,
    # store this column name
    target_col_name <- target_names
  } else {
    target_col_name <- NULL
  }
  if (!is.null(df)) {
    get_args_from_df(df, args = c("target_names", "F", "Q", "w", "kappa", "alpha", "g", "dg"))
    N <- nrow(df)
  } else {
    N <- length(F)
  }
  # if target_names still not provided, name targets with numbers
  if ((length(target_names) == 1) && (is.na(target_names))) {
    target_names <- as.character(1:N)
  }
  # name the target weights
  if (length(w) == 1) {
    w <- rep(w, N) %>% set_names(target_names)
  } else if (length(w) == N) {
    w <- w %>% set_names(target_names)
  } else {
    stop("w wrong length")
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
  # Old initiation code that doesn't seem useful: the algorithm will go to medians
  # on fist step whether Infs are there or not.
  # qs[qs==Inf] <- (w^(-1)*rep(1,N)*2*max(K))[qs==Inf]

  # get marginal expected benefits
  Lambda <- meb_gpl_df(df = gpl, F = F, w = w)
  # store their values at boundary
  Lambda0 <- map2_dbl(Lambda, 0, exec)
  # take max as an upper bound on lambda
  lamU0 <- max(Lambda0)
  # set up a data frame to track iterations for all K
  # with columns for binary search interval endpoints at lambda = MEB = 0 (corresponding to
  # allocation at quantile) and the maximum over target set of MEB's evaluated at x = 0
  # (corresponding to the maximum MEB if the allocations were all reduced as much as possible - to 0)
  # as well as columns for the current lambda and previous lambda iterates
  Kdf <- tibble(
    K = K,
    xs = list(gpl %>% select(tidyselect::all_of(target_col_name)) %>% mutate(`0` = qs)), # for iterations
    x = list(tibble::deframe(xs[[1]])), # for final allocation
    qs_OK = map(x, ~sum(w*.)) < K,
    converged = FALSE,
    lamL = 0,
    lamL_seq = list(0),
    lamU = lamU0,
    lamU_seq = list(lamU0),
    lam = 0,
    lam_seq = list(0),
    lam_prev = 0
  ) %>% dplyr::arrange(K)
  # remove rows where quantiles satisfy constraint
  out <- Kdf %>% dplyr::filter(qs_OK) %>% mutate(converged = TRUE)
  Kdf <- Kdf %>% dplyr::filter(!qs_OK)
  # return x = 0 if there is no marginal benefit to allocating more than 0 in any component
  if ((nrow(Kdf) > 0 ) && (Kdf$lamU[1] <= 0)) {
    Kdf$x <- list(rep(0, N))
    out <- rbind(out, Kdf) %>% arrange(K) # would out ever be non-empty here?
    message("All targets receiving zero allocation")
    return(out)
  }
  # create counter for labeling iterates
  tau <- 1
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
      K_lam <- max(Kdf_lam$K)
      x_tau <- Kdf_lam$x[[1]] # this is the whole point of using Kdf
      stopifnot(length(lam_prev <- unique(Kdf_lam$lam_prev)) == 1) # sanity check
      for (i in 1:N) {
        # if (tau == 13 && i == 1) browser()
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
            if (is.finite(x_tau[i])) {
              # Use x_tau to define search intervals; motivation here is
              # to make it harder for errors to occur silently
              if (lam_tau < lam_prev) {
                # if lam has decreased we need to look closer to the quantile
                x_tau[i] <- unirootL(
                  f = lam_grad,
                  lower = x_tau[i],
                  upper = K_lam,
                  point_mass_window = point_mass_window
                )
              } else {
                # and if not we need to look further toward 0
                x_tau[i] <- unirootL(
                  f = lam_grad,
                  lower = 0,
                  upper = x_tau[i] + point_mass_window,
                  point_mass_window = point_mass_window
                )
              }
            } else {
              # if x_tau = Inf (as for alpha = 1) we need to use expanding intervals
              if (lam_tau < lam_prev) {
                stop("lambda should not decrease while there are infinite allocations")
              } else {
                x_tau[i] <- uniroot(
                  f = lam_grad,
                  lower = 0,
                  upper = K_lam * (1 + point_mass_window),
                  extendInt = "downX"
                )$root
              }
            }
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
    x = pmap(list(post_processed, x, K, lam),
             function(post_processed, x, K, lam) {
               if (post_processed) {
                 x <- post_process(x, K, lam, w, Lambda, eps_lam, point_mass_window)
               }
               return(x)
             }
             ))
  # organize
  out <- out %>% select(-c(converged, lam_prev)) %>% arrange(K)
  # create a scorable data frame for each K
  out <- out %>% mutate(xdf = map(x, function(x) {
    tibble::enframe(x, name = target_col_name, value = "x") %>%
      dplyr::left_join(gpl[, c(target_col_name, "gpl_loss_fun")], by = target_col_name) %>%
      mutate(score_fun = map2(x, gpl_loss_fun, function(x, gpl_loss_fun) {
        function(y) {gpl_loss_fun(x, y)}
      })) %>% select(-gpl_loss_fun)
  }))
  return(structure(
    out,
    class = c("allocated", class(out)),
    gpl_df = gpl,
    w = w,
    target_col_name = target_col_name))
}

#' Post process allocations that are infeasible due to plateaus
#'
#' @param x last iteration of allocation via lambda search
#' @param K constraint
#' @param lam last lambda iterate
#' @param w weights
#' @param Lambda meb's
#' @param eps_lam lambda tolerance that led to lambda convergence failure
#' @param point_mass_window
#'
#' @return
#' @export
#'
#' @examples
post_process <- function(x, K, lam, w, Lambda, eps_lam, point_mass_window) {
  if (sum(w * x) > K) {
    lam_grad_eps <- map(Lambda, ~ function(x) .(x) - lam - eps_lam)
    x_L <- x_U <- x
    for (i in 1:length(x)) {
      if (lam_grad_eps[[i]](0) <= 0 || x[i] < point_mass_window) {
        x_L[i] <- 0
      } else {
        tryCatch({x_L[i] <- unirootL(
          f = lam_grad_eps[[i]],
          lower = 0,
          upper = x[i],
          point_mass_window = point_mass_window
        )}, error = function(e) {
          message(paste(
            "Error at index =", i,
            ", K =", K,
            ": ", e$message))
          return(NULL)
        } )
      }
    }
  } else {
    # We want to catch approximate right endpoints of plateaus but this will
    # not be possible when lam is close to zero (because K large) - i.e., when
    # we are on the upper tails of the forcasts.
    # In that case we just look for x_U's larger than x_L's (which finding a root
    # for lam/2 should accomplish) and let uniroot extend
    # the t search interval to the right in the final step below involving
    # Delta.
    lam_grad_eps <- map(Lambda, ~ function(x) .(x) - max(lam - eps_lam, lam/2))
    x_L <- x_U <- x
    for (i in 1:length(x)) {
      if (lam_grad_eps[[i]](0) <= 0) {
        x_U[i] <- x[i]
      } else {
        tryCatch({x_U[i] <- uniroot(
          f = lam_grad_eps[[i]],
          lower = x[i],
          upper = x[i] + 1,
          extendInt = "downX"
        )$root}, error = function(e) {
          message(paste(
            "Error at index =", i,
            ", K =", K,
            ": ", e$message))
          return(NULL)
        } )
      }
    }
  }
  Delta <- function(t) {
    sum(w * ((1 - t) * x_L + t * x_U)) - K
  }
  if ((Delta(0) < 0) & (Delta(1) > 0)) {
    t_star <- uniroot(f = Delta, interval = c(0, 1))$root
    return((1 - t_star) * x_L + t_star * x_U)
  } else if (lam <= 2 * eps_lam) {
    t_star <- uniroot(f = Delta, interval = c(0, 1), extendInt = "upX")$root
    return((1 - t_star) * x_L + t_star * x_U)
  } else {
    stop("Post-processing failed")
  }
}

#' Get the `gpl_df` attribute of an allocated data frame which containes gpl loss functions
#' for each target and the parameters used to form them.
#'
#' @param adf an allocated data frame
#'
#' @return data frame
#' @export
#'
#' @examples
gpl <- function(adf) {
  if (!inherits(adf, "allocated")) {
    stop("Input must be of class 'allocated'")
  }
  attr(adf, "gpl_df")
}

#' Get the weights of an allocated data frame
#'
#' @param adf an allocated data frame
#'
#' @return weights used in allocation
#' @export
#'
#' @examples
weights.allocated <- function(adf) {
  attr(adf, "w")
}

#' Make an allocated data frame slim, i.e., containing
#' columns only for scores and optionally a list column or unnested
#' columns for the `xdf` component data frames of the allocation
#'
#' @param adf an allocated data frame
#' @param xdf_action character indicating whether to unnest `xdf`.
#' Default is to unnest when data from is scored and leave nested
#' when data frame is scored.
#' @param id_cols columns to also keep such as model name and origin time
#' @param rm_score_fun_if_scored remove the scoring function in the xdf list column of a scored
#'  data frames to save memory; defaults to TRUE since a scoring function is usually not needed for a
#'  scored data frame.
#' @param rm_score_fun_if_not_scored defaults to FALSE since a scoring function can be used to
#' score the data frame.
#'
#' @return a slim allocated data frame
#' @export
#'
#' @examples
slim <- function(
    adf,
    xdf_action = c("default", "unnest", "nest"),
    id_cols = NULL,
    rm_score_fun_if_scored = TRUE,
    rm_score_fun_if_not_scored = FALSE) {
  # stopifnot("allocated" %in% class(adf))
  if (!is.null(id_cols)) {
    adf <- adf %>% dplyr::relocate(id_cols)
  }
  id_cols <- c(id_cols, "K", "xdf")
  xdf_action <- match.arg(xdf_action)
  if ("scored" %in% class(adf)) {
    out <- adf %>%
      select(tidyselect::any_of(c(
        id_cols, "score_raw", "score_oracle", "score"
      )))
    if (rm_score_fun_if_scored) {
      out <- out %>% mutate(xdf = map(xdf, ~select(.,-tidyselect::any_of("score_fun"))))
    }
    if (xdf_action == "unnest") {
      out <- out %>% tidyr::unnest(xdf)
    }
  } else {
    out <- adf %>% select(tidyselect::any_of(id_cols))
    if (rm_score_fun_if_not_scored) {
      out <- out %>% mutate(xdf = map(xdf, ~select(.,-tidyselect::any_of("score_fun"))))
    }
    if (xdf_action != "nest") {
      out <- out %>% tidyr::unnest(xdf)
    }
  }
  class(out) = c("slim", class(adf))
  return(structure(
    out,
    gpl_df = attr(adf, "gpl_df"),
    w = attr(adf, "w"),
    target_col_name = attr(adf, "target_col_name")
  ))
}

#' Score the allocations from set of forecasts against realized outcomes `y`.
#'
#' @param df data frame containing columns for other list arguments which are used only when
#'  those arguments (e.g. `F`) are not passed directly
#'
#' @return a data frame of the form returned by `allocate` with additional columns for
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
#' @param slim logical; if TRUE use slim data frames
#' @inheritParams allocate
#' @rdname alloscore
#' @export
alloscore.default <- function(df = NULL, K, target_names = NA,
                              y, F = NULL, Q = NULL, w = 1,
                              kappa = 1,
                              alpha = 1,
                              g = "x",
                              dg = NA,
                              eps_K = .01,
                              eps_lam = 1e-4,
                              against_oracle = TRUE,
                              slim = FALSE) {

  # allocate will handle validation and attribute assignments
  # if (!is.null(df)) {
  #   get_args_from_df(df)
  # }
  a <- allocate(
    df = df,
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
  )
  y <- setNames(y, nm = gpl(a)[[attr(a, "target_col_name")]])
  if (slim) {a <- slim(a)}

  return(alloscore(a, y = y, against_oracle = against_oracle))
  # uses allocated or slim method
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
  scored_df <- df %>% mutate(
    xdf = map(xdf, function(.xdf) {
      .xdf %>% mutate(
        y = y,
        components_raw = map2_dbl(score_fun, y, function(.score_fun, .y) {.score_fun(.y)})
        )
    }),
    ytot = map_dbl(xdf, ~sum(attr(df, "w") * dplyr::pull(., y))),
    score_raw = map_dbl(xdf, ~sum(dplyr::pull(., components_raw)))
  ) %>% relocate(score_raw, ytot, .after = K)

  if (against_oracle) {
    # move this code to oracle_alloscore?
    oracle_scores <- gpl(df) %>%
      oracle_allocate(y = y,
                      w = attr(df, "w"),
                      K = df$K,
                      target_names = attr(df, "target_col_name")) %>%
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
      ) %>% relocate(score, .after = K) %>% relocate(score_oracle, .after = score_raw)
  }
  class(scored_df) <- c("scored", class(scored_df))
  return(scored_df)
}

#' Score a slim allocated data frame
#'
#' @param slim_df a slim allocated data frame
#' @param ys a list of named outcome vectors with names matching target names in slim_df
#'
#' @return a tibble with columns `xdf` and `scores` containing the components and scores for
#'  each sample in `ys`
#' @export
#'
#' @examples
alloscore.slim <- function(slim_df, ys, against_oracle) {
  stopifnot("allocated" %in% class(slim_df))
  if (!is.list(ys)) {ys <- list(ys)}
  if (!all(purrr::map_lgl(ys,
                          ~ !is.null(names(.)) &
                          identical(names(ys[[1]]), names(.))
                          ))) {
    stop("ys need to be consistently named by their targets")
  }
  target_col_name <- attr(slim_df, "target_col_name")
  comp_base <- slim_df %>% select(K, target_col_name, x)
  w <- weights(slim_df)
  Ks <- unique(slim_df$K)
  gpl_fns <- gpl(slim_df)$gpl_loss_fun
  scores <- map(1:length(ys), function(i) {
    y <- ys[[i]]
    oracle_cols <- map_df(Ks, ~oracle_alloscore_direct(y, ., w, gpl_fns))
    xdf <- bind_cols(comp_base, y = rep(y, length(Ks)), oracle_cols)
    xdf$components_raw <- map2_dbl(
      slim_df[[target_col_name]],
      slim_df$score_fun,
      function(nm, fn) {fn(y[nm])})
    xdf$components <- xdf$components_raw - xdf$components_oracle
    scores <- xdf %>% group_by(K) %>%
      summarise(
        ytot = sum(w * y),
        score_raw = sum(components_raw),
        score_oracle = sum(components_oracle),
        score = sum(components))
    return(list(xdf = xdf, scores = scores))
  }, .progress = TRUE)
  out <- tibble::as_tibble(purrr::transpose(scores)) %>%
    rownames_to_column() %>% rename(samp = rowname)
  return(out)
}

#' Allocate according to an oracle's knowledge of outcome y
#'
#' @param gpl_df a `gpl_df` object
#' @param y observed outcomes
#' @inheritParams allocate
#'
#' @return see `allocate`
#' @export
#'
#' @examples
oracle_allocate <- function(gpl_df, y, K, w, ...) {
  if (inherits(gpl_df, "allocated")) {
    w <- weights(gpl_df)
    gpl_df <- gpl(gpl_df)
  }
  # not dealing with target_name_col identification for now
  gpl_df <- gpl_df %>% select(g:alpha) %>%
    mutate(kappa = alpha, alpha = 1)
  return(allocate(
    df = gpl_df,
    w = w,
    F = map(y, function(y) {function(x) {1*(x>=y)}}),
    Q = map(y, function(y) {function(p) {y}}),
    K = K,
    eps_lam = .01, # prevent likely unnecessary search steps
    ...))
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
                             g = "x", ...) {
  stop("Not implemented yet")
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

#' Directly find the oracle allocation with no data management.
#' Assumes for now that all targets have same gpl parameters and
#' (probably?) that `g = x`
#'
#' @param y
#' @param K
#' @param w
#'
#' @return
#' @export
#'
#' @examples
oracle_allocate_direct <- function(y, K, w) {
  if (sum(w * y) <= K) {
    return(y)
  } else {
    return(y * K / sum(w * y))
  }
}

#' Directly find and score the oracle allocation with no data management.
#'
#' @param y
#' @param K
#' @param w
#' @param gpl_fns
#'
#' @return
#' @export
#'
#' @examples
oracle_alloscore_direct <- function(y, K, w, gpl_fns) {
  oracle <- oracle_allocate_direct(y, K, w)
  components_oracle <- pmap_dbl(
    list(gpl_fns, oracle, y), function(fn, o, y) {fn(o, y)})
  return(tibble(oracle = oracle, components_oracle = components_oracle))
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
