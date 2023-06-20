

#' Plot components generic
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_components <- function(...) {
  UseMethod("plot_components")
}

#' Plot components of a nested data frame of scored allocations
#'
#' @param df
#' @param Ks
#' @param scored_col_name
#' @param origin_time_col_name
#' @param model_col_name
#' @param target_col_name
#'
#' @return
#' @export
#'
#' @examples
plot_components.default <- function(
    df,
    Ks = NULL,
    scored_col_name = "scored",
    origin_time_col_name = NULL,
    model_col_name = NULL,
    target_col_name = NULL) {
  scored_list <- df[[scored_col_name]]
  if (!is.null(origin_time_col_name)) {
    origin_time <- df[[origin_time_col_name]]
    scored_list <- map2(
      origin_time, scored_list,
      function(.origin_time, .scored_list) {
        .scored_list[[origin_time_col_name]] <- .origin_time
        return(.scored_list)
      })
  }
  if (!is.null(model_col_name)) {
    model <- df[[model_col_name]]
    scored_list <- map2(
      model, scored_list,
      function(.model, .scored_list) {
        .scored_list[[model_col_name]] <- .model
        return(.scored_list)
      })
  }
  # scored_list <- purrr::pmap(df,
  #                           function(
  #   origin_time_col_name, model_col_name, scored_col_name, ...) {
  #                             alloc %>% mutate(origin_time = origin_time,
  #                                              model = model,
  #                                              .before = 1)
  #                           })
  do.call(plot_components, c(
    scored_list,
    Ks = Ks,
    origin_time_col_name = origin_time_col_name,
    model_col_name = model_col_name,
    target_col_name = target_col_name
    ))
}

#' Default plot_components method; used to plot multiple allocated dfs.
#'
#' @param ... allocated data frames
#' @param Ks Ks to plot over, defaults to NULL which keeps all. Note that some plots
#'  require a single K
#' @param origin_time_col_name name of of origin time columns
#' @param model_col_name name of model columns
#' @param target_col_name name of columns containing target names
#'
#' @return
#' @export
#'
#' @examples
plot_components.allocated <- function(
    ...,
    Ks = NULL,
    origin_time_col_name = NULL,
    model_col_name = NULL,
    target_col_name = NULL) {
  adfs <- list(...)
  reqs <- map(
    list(
      is_allocated = ~ "allocated" %in% class(.),
      is_scored = ~ "scored" %in% class(.),
      has_origin_time = ~ !is.null(origin_time_col_name) && origin_time_col_name %in% names(.),
      has_model = ~ !is.null(model_col_name) && model_col_name %in% names(.)
    ),
    ~purrr::map_lgl(adfs, .)
  )

  if (!all(reqs$is_allocated)) {
    stop("Data frames must be allocated.")
  }
  if (!all(reqs$is_scored)) {
    stop("Data frames must be scored.")
  }
  if (any(reqs$has_origin_time) && any(!reqs$has_origin_time)) {
    stop("Data frames must have a consistent origin time column.")
  }
  if (any(reqs$has_model) && any(!reqs$has_model)) {
    stop("Data frames must have a consistent model column.")
  }
  if (!any(reqs$has_model) && !any(reqs$has_origin_time)) {
    # assume that each df for a different model
    for (i in 1:length(adfs)) {
      adfs[[i]]$model <- paste0("model_", i)
    }
    model_col_name <- "model"
  }
  adf <- bind_rows(adfs)
  key <- adf %>% select(tidyr::any_of(c(origin_time_col_name, model_col_name, "K")))
  if (any(duplicated(key))) {
    stop("Allocations are not uniquely identified.")
  }
  slim_df <- slim(adf, id_cols = c(origin_time_col_name, model_col_name))
  return(plot_components.slim(
    slim_df, Ks = Ks,
    origin_time_col_name = origin_time_col_name,
    model_col_name = model_col_name,
    target_col_name = target_col_name))
}

#' Plot components
#'
#' @param slim_df
#' @param Ks
#' @param origin_time_col_name
#' @param model_col_name
#' @param target_col_name
#'
#' @return
#' @export
#'
#' @examples
plot_components.slim <- function(
    slim_df,
    Ks = NULL,
    origin_time_col_name = NULL,
    model_col_name = NULL,
    target_col_name = NULL) {
  if (!"scored" %in% class(slim_df)) {
    stop("data frames must be scored")
  }
  if (!is.null(Ks)) {
    slim_df <- slim_df %>% filter(K %in% Ks)
  }
  if ("xdf" %in% names(slim_df)) {
    slim_df <- slim_df %>% tidyr::unnest(xdf)
  }

  if (!is.null(origin_time_col_name) && origin_time_col_name %in% names(slim_df)) {
    slim_df <- slim_df %>% rename(origin_time = !!rlang::sym(origin_time_col_name))
  }

  if (!is.null(model_col_name) && model_col_name %in% names(slim_df)) {
    slim_df <- slim_df %>% rename(model = !!rlang::sym(model_col_name))
  }

  if (!"model" %in% names(slim_df)) {
    slim_df$model <- "model"
  }

  if (!is.null(target_col_name) && target_col_name %in% names(slim_df)) {
    slim_df <- slim_df %>% rename(target_names = !!rlang::sym(target_col_name))
  }
  slim_df <- slim_df %>% bind_rows(
  slim_oracle <- slim_df %>% filter(model == model[1]) %>%
    mutate(model = "oracle",
           components_raw = components_oracle)
  )
  if (length(unique(slim_df$K)) == 1) {
    p <- ggplot(slim_df, aes(x = model,
                             y = components_raw,
                             fill = target_names)) +
      geom_bar(position = "stack", stat = "identity")
    if ("origin_time" %in% names(slim_df)) {
      p <- p + facet_wrap( ~ origin_time)
    }
  } else {
    p <- ggplot(slim_df, aes(x = K,
                             y = components_raw,
                             fill = target_names)) +
      geom_area()
    if ("origin_time" %in% names(slim_df) && "model" %in% names(slim_df)) {
      p <- p + facet_grid(rows = vars(model), cols = vars(origin_time))
    }
    if ("origin_time" %in% names(slim_df) && !"model" %in% names(slim_df)) {
      p <- p + facet_grid(cols = vars(origin_time))
    }
    if (!"origin_time" %in% names(slim_df) && "model" %in% names(slim_df)) {
      p <- p + facet_grid(rows = vars(model))
    }
  }
  return(p)
}


