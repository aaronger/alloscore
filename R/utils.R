#' Produce a pdqr function with set parameters
#'
#' @param ... list of parameter name-value pairs to be partialed out
#' @param dist distribution root name or `distfromq`
#' @param type p, d, q, or r
#' @inheritParams add_pdqr_funs
#'
#' @return function of type p, d, q, or r
#' @export
#'
#' @examples
pdqr_factory <- function(..., dist, type,
                         trans = NULL, trans_inv = NULL,
                         transpars = NULL) {
  pars <- list(...) %>% .[!is.na(.)]
  if (dist == "distfromq") {
    funfac <- get(paste0("make_", type, "_fn"))
    selected_pars <- pars[intersect(names(formals(funfac)), names(pars))]
    Fun <- exec(funfac, !!!selected_pars)
  } else {
    statfun <- get(paste0(type, dist))
    selected_pars <- pars[intersect(names(formals(statfun)), names(pars))]
    Fun <- partial(statfun, ...=, !!!selected_pars)
  }
  if (!is.null(trans)) {
    transpars <- c(pars[intersect(names(formals(trans)), names(pars))], transpars)
    if (type %in% c("q", "r")) {
      Fun_comp <- compose(partial(trans_inv, ...=, !!!transpars), Fun)
    } else {
      Fun_comp <- compose(Fun, partial(trans, ...=, !!!transpars))
    }
    Fun <- Fun_comp
  }
  return(Fun)
}

#' Add columns for pdqr functions (i.e., functions to evaluate the cdf, pdf, or
#' quantile function, or to generate random deviates) determined by parameter columns
#'
#' @param df data frame; must have columns for each parameter required by any
#'  of the distribution functions indicated in a column named `dist` provided directly
#'  or via argument `dist`
#' @param dist character vector of length 1 or `nrow(df)` containing distribution
#'  root names (or `distfromq`); must be provided if df does not have a column `dist`
#' @param types which pdqr function list columns to add
#' @param trans function with arguments x for outcome and parameters for a transformation
#'  precomposed with pd functions.
#' @param trans_inv inverse of trans with argument p for probability level; will be
#'  post-composed with qr function.
#' @param transpars parameters for trans and trans_inv not included in columns of df
#' @param fnames names of pdqr function list columns
#'
#' @return a tibble with list columns containing selected pdqr functions
#' @export
#'
#' @examples
add_pdqr_funs <- function(
    df,
    dist = df$dist,
    trans = NULL,
    trans_inv = NULL,
    transpars = NULL,
    types = c("p", "d", "q", "r"),
    fnames = c(p = "F", d = "f", q = "Q", r = "r")) {
  if (!"dist" %in% names(df)) {
    if (is.null(dist)) stop("distributions must be specified")
    df$dist <- dist
  }
  if (!"trans" %in% names(df) & !is.null(trans)) {
    if (!is.list(trans)) {trans <- list(trans)}
    df$trans <- trans
  }
  if (!"trans_inv" %in% names(df) & !is.null(trans_inv)) {
    if (is.null(trans_inv) & !is.null(trans) & any(is.element(c("q", "r"), types))) {
      stop("inverse transformation required for q or r functions")
    }
    if (!is.list(trans_inv)) {trans_inv <- list(trans_inv)}
    df$trans_inv <- trans_inv
  }
  for (type in types) {
    df <- df %>% mutate(
      !!fnames[type] := pmap(., .f = pdqr_factory, type = type)
    )
    }
  return(as_tibble(df))
}

#' Convert newsvendor parameters to kappa, alpha
#'
#' @param ax wholesale cost
#' @param ay
#' @param a_minus
#' @param a_plus
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_news_params <- function(ax, ay = NULL, a_minus, a_plus) {
  return(data.frame(
    kappa = as.double(a_plus + a_minus),
    alpha = (a_minus - ax) / (a_plus + a_minus)
  ))
}

#' Convert over/underprediction parameters to kappa, alpha
#'
#' @param O incremental loss that is incurred when overprediction leads to unused supply
#' @param U incremental loss that is incurred when underprediction leads to unmet need
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_ou_params <- function(O, U) {
  return(data.frame(
    kappa = O + U,
    alpha = U/(O+U)
  ))
}

#' Convert meteorologist parameters to kappa, alpha
#'
#' @param C marginal cost per unit of recommended protection
#' @param L marginal loss due to under-provision of needed resources
#'
#' @return kappa and alpha as list
#' @export
#'
#' @examples
stdize_met_params <- function(C, L) {
  return(data.frame(
    kappa = L,
    alpha = 1-C/L
  ))
}

#' Utility function to make data frames of forecasts easier to work with;
#' works via side-effects of `assign`.
#' Note that we need the "empty" arguments to actually have NULL or NA defaults since
#' leaving them actually empty seems to create un-inspectable promises when the call stack
#' goes beyond just the caller of this function.
#' @param df data frame with columns we might want caller to assign to its arguments
#' @param args character list of arguments we want df to potentially provide
get_args_from_df <- function(df, args) {
  e <- rlang::caller_env()
  # get arguments that have non-NULL defaults
  has_default_in_caller <- names(purrr::discard(formals(rlang::caller_fn()), is.null))
  has_default <- intersect(has_default_in_caller, args)
  # overwrite if they are supplied by df
  replacements <- intersect(has_default, names(df))
  for (name in replacements) {
    assign(name, df[[name]], envir = e)
  }
  # supply arguments from df when the default is NULL
  empty_objects <- purrr::keep(as.list(e), ~(length(.) < 2 && (is.null(.) || is.na(.))))
  empty_arg_names <- intersect(names(empty_objects), args)
  supplied_by_df <- intersect(empty_arg_names, names(df))
  for (name in supplied_by_df) {
    assign(name, df[[name]], envir = e)
  }
  # if value of arg is a string specifying a column of df,
  # assign that column to arg
  string_args <- c()
  for (name in setdiff(args, empty_arg_names)) {
    arg_val <- e[[name]]
    if (is.character(arg_val) && (length(arg_val) == 1) && (arg_val %in% names(df))) {
      assign(name, df[[arg_val]], envir = e)
      string_args <- c(name, string_args)
    }
  }
  # make apparent to caller what came from df
  assign("args_from_df", c(replacements, supplied_by_df, string_args), envir = e)
}

#' Get a function of x from a string, e.g., "log(x)"
get_function <- function(func) {
  # Ensure 'x' is used as the variable name
  if (!grepl("\\bx\\b", func)) {
    stop("Variable 'x' must be used in the function string")
  }
  if (func == "x") {
    return(function(x) x)
  }
  if (func == "log(x)") {
    return(function(x) log(x))
  }
  if (func == "exp(x)") {
    return(function(x) exp(x))
  }
  func_expr <- parse(text = func)
  return(function(x) eval(func_expr))
}

#' Get the derivative of a function of x entered as a string, e.g., "log(x)"
get_derivative <- function(func) {
  # Ensure 'x' is used as the variable name
  if (!grepl("\\bx\\b", func)) {
    stop("Variable 'x' must be used in the function string")
  }
  if (func == "x") {
    return(function(x) 1)
  }
  if (func == "log(x)") {
    return(function(x) 1/x)
  }
  if (func == "exp(x)") {
    return(function(x) exp(x))
  }
  # Turn string into an expression
  expr <- parse(text = func)
  # Differentiate
  der_expr <- D(expr, "x")
  # Convert to function
  der_func <- function(x) eval(der_expr, envir = list(x = x))
  return(der_func)
}
