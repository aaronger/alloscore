

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
    target_col_name = NULL,
    show_oracle = TRUE) {
  scored_list <- df[[scored_col_name]]
  if (is.null(scored_list)) {
    stop("No scored data frames were specified. If given in a column, pass column name  to 'scored_col_name'")
  }
  # attach origin_time and model columns to each scored data frame
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
    target_col_name = target_col_name,
    show_oracle = show_oracle
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
    target_col_name = NULL,
    show_oracle = TRUE) {
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
    if (length(adfs) > 1) {
      stop("Please include an origin time or model column in allocated data frames.")
    }
    # for (i in 1:length(adfs)) {
    #   adfs[[i]]$model <- paste0("model_", i)
    # }
    # model_col_name <- "model"
  }
  adf <- bind_rows(adfs)
  key <- adf %>% select(tidyr::any_of(c(origin_time_col_name, model_col_name, "K")))
  if (any(duplicated(key))) {
    stop("Allocations are not uniquely identified.")
  }
  slim_df <- slim(adf, id_cols = c(origin_time_col_name, model_col_name))
  return(plot_components_slim(
    slim_df, Ks = Ks,
    origin_time_col_name = origin_time_col_name,
    model_col_name = model_col_name,
    target_col_name = target_col_name,
    show_oracle = show_oracle))
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
plot_components_slim <- function(
    slim_df,
    Ks = NULL,
    origin_time_col_name = NULL,
    model_col_name = NULL,
    target_col_name = NULL,
    show_oracle = TRUE,
    show_raw = TRUE,
    order_at_K = NULL,
    order_at_model = NULL,
    order_at_origin_time = NULL,
    pal_top = c("#CD3333", "#009ACD", "#FFB90F", "#8B3E2F", "#8B8B00"),
    bar_positioning = "stack"
    ) {
  # if (!"scored" %in% class(slim_df)) {
  #   stop("data frames must be scored")
  # }
  if (!is.null(Ks)) {
    slim_df <- slim_df %>% filter(K %in% Ks)
  }
  if ("xdf" %in% names(slim_df)) {
    slim_df <- slim_df %>% tidyr::unnest(xdf)
  }
  # normalize column names
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
  # add on an oracle data frame
  if (show_oracle) {
    slim_df <- slim_df %>% bind_rows(
    slim_oracle <- slim_df %>% filter(model == model[1]) %>%
      mutate(model = "oracle",
             components_raw = components_oracle)
    )
  }
  if (show_raw) {
    slim_df <- slim_df %>% mutate(components = components_raw)
  }
  # order the components as factors
  lev_df <- slim_df
  if (!is.null(order_at_K)) {
    lev_df <- lev_df %>% filter(K == order_at_K)
  }
  if (!is.null(order_at_model)) {
    lev_df <- lev_df %>% filter(model == order_at_model)
  }
  if (!is.null(order_at_origin_time)) {
    lev_df <- lev_df %>% filter(K == order_at_origin_time)
  }
  levs <- lev_df %>% mutate(
    target_names = forcats::fct_reorder(target_names, desc(components))) %>%
    pull(target_names) %>% levels()

  if (!is.null(order_at_K) || !is.null(order_at_model) || !is.null(order_at_origin_time)) {
    colors <- setNames(c(pal_top, rep(c("lightgrey", "darkgrey"), length.out = 60)), levs)
  } else {
    colors <- setNames(palette.colors(n = 60, palette = "Paired", recycle = T), levs)
  }
  slim_df <- slim_df %>% mutate(target_names = factor(target_names, levels = levs))

  ### Main plotting code
  # plot bars and facet by time and model if there's only 1 K
  if (length(unique(slim_df$K)) == 1) {
    if (bar_positioning == "stack") {
    p <- ggplot(slim_df, aes(x = model,
                             y = components,
                             fill = target_names)) +
      geom_bar(position = "stack", stat = "identity")
    } else if (bar_positioning == "dodge") {
      p <- ggplot(slim_df, aes(x = target_names,
                               y = components,
                               fill = model)) +
        geom_bar(position = "dodge", stat = "identity")
    }
    if ("origin_time" %in% names(slim_df)) {
      p <- p + facet_wrap( ~ origin_time)
    }
  # plot ribbons if there are multiple Ks
  } else {
    p <- ggplot(slim_df, aes(x = K,
                             y =  components,
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
  if (bar_positioning == "stack") {
    p <- p + scale_fill_manual(values = colors)
  }
  return(p)
}


#' Plot allocation scores
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
plot_scores_slim <- function(
    slim_df,
    Ks = NULL,
    ts = FALSE,
    ts_dates = NULL,
    origin_time_col_name = "reference_date",
    model_col_name = "model",
    target_col_name = "abbreviations",
    show_oracle = FALSE,
    show_raw = FALSE,
    order_at_K = NULL,
    order_at_model = NULL,
    order_at_origin_time = NULL,
    palette = score_palette,
    linetypes = score_linetypes,
    ytot = TRUE
  ) {
  # if (!"scored" %in% class(slim_df)) {
  #   stop("data frames must be scored")
  # }
  if (!is.null(Ks)) {
    slim_df <- slim_df %>% filter(K %in% Ks)
  }
  if ("xdf" %in% names(slim_df)) {
    slim_df <- slim_df %>% tidyr::unnest(xdf)
  }
  # normalize column names
  slim_df <- normalize_col_names(
    slim_df,
    origin_time_col_name,
    model_col_name,
    target_col_name)
  # add on an oracle data frame
  if (show_oracle) {
    slim_df <- slim_df %>% bind_rows(
    slim_oracle <- slim_df %>% filter(model == model[1]) %>%
      mutate(model = "oracle",
             components_raw = components_oracle)
    )
  }

  # order the models by score
  lev_df <- slim_df
  if (!is.null(order_at_K)) {
    lev_df <- lev_df %>% filter(K == order_at_K)
  }
  if (!is.null(order_at_model)) {
    lev_df <- lev_df %>% filter(model == order_at_model)
  }
  if (!is.null(order_at_origin_time)) {
    lev_df <- lev_df %>% filter(K == order_at_origin_time)
  }
  levs <- lev_df %>% mutate(
    model = forcats::fct_reorder(model, desc(score))) %>%
    pull(model) %>% levels()

  slim_df <- slim_df %>% mutate(model = factor(model, levels = levs))

  ### Main plotting code
  # plot against time time and model if there's only 1 K
  if ((length(unique(slim_df$K)) == 1) ||
      (ts && length(Ks) == length(unique(slim_df$origin_time)))) {
    ts_base <- ts_base <- tidyr::expand_grid(
      model = mkeep,
      tibble(origin_time = sort(unique(slim_df$origin_time)),
             K = Ks))
    slim_df_ts <- ts_base %>% dplyr::left_join(slim_df,
                        by = c("model", "origin_time", "K")) %>%
      group_by(origin_time, model) %>% slice(1) %>% ungroup()
    if (!is.null(ts_dates)) {
      slim_df_ts <- slim_df_ts %>% filter(origin_time %in% ts_dates)
    }
    slim_df_ts <- slim_df_ts %>% mutate(origin_time = as.Date(origin_time))
    p <- ggplot(slim_df_ts, aes(x = origin_time,
                             y = score,
                             group = model,
                             color = model,
                             linetype = model)) +
      geom_line() + geom_point()
  # plot against K
  } else {
    p <- ggplot(slim_df, aes(x = K,
                             y =  score,
                             color = model,
                             linetype = model)) +
      geom_line()
    if ("origin_time" %in% names(slim_df) && "model" %in% names(slim_df)
        && unique(slim_df[["origin_time"]] > 1 )) {
      p <- p + facet_grid(cols = vars(origin_time))
    }
    if ("origin_time" %in% names(slim_df) && !"model" %in% names(slim_df)) {
      p <- p + facet_grid(cols = vars(origin_time))
    }
    if (!"origin_time" %in% names(slim_df) && "model" %in% names(slim_df)) {
      p <- p + facet_grid(rows = vars(model))
    }
  }
  p <- p + scale_color_manual(values = palette) +
    scale_linetype_manual(values = linetypes)

  if (ytot) {
    ytot_df <- slim_df %>% group_by(origin_time, model, K) %>%
        summarise(ytot = sum(y))
    p <- p + geom_vline(
      data = ytot_df,
      aes(xintercept = ytot),
      linetype = 2)
  }
  return(p)
}

#' Column name normalizing utility
#'
#' @param slim_df
#' @param origin_time_col_name
#' @param model_col_name
#' @param target_col_name
#'
#' @return
#'
#' @examples
normalize_col_names <- function(
    slim_df,
    origin_time_col_name,
    model_col_name,
    target_col_name
    ) {
  # normalize column names
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
  slim_df
}

