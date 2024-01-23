

#' Title Additional code to generalize the results
#'
#' @param SL.library the learners we input
#' @param n_split number of splits
#'
#' @return it return more general results
#'
#' @export
#'
#' @examples
#'
#' # See the README file
#'
generate_placeholder_output <- function(SL.library, n_split) {
  placeholder <- list()

  placeholder$results <- tibble::tibble(
    rd = NA,
    var = NA
  )

  placeholder_weights <- tibble::tibble(row = 1:(2 * n_split))

  for (item in SL.library) {
    item_col <- paste0(item, "_All")
    placeholder_weights <- dplyr::mutate(placeholder_weights, !!item_col := 0)
  }

  model_values <- rep(c("x", "y"), each = n_split)
  split_values <- rep(1:n_split, 2)
  prev_values <- rep(rep(NA, 2 * n_split), length.out = 2 * n_split)

  placeholder_weights <- dplyr::mutate(
    placeholder_weights,
    model = factor(model_values),
    split = split_values,
    prev = prev_values
  )

  placeholder_weights <- dplyr::select(placeholder_weights, -row)

  placeholder$weight <- placeholder_weights

  return(placeholder)
}


summarize_multiple <- list(
  mean = ~mean(.x, na.rm = TRUE),
  max = ~max(.x, na.rm = TRUE),
  min = ~min(.x, na.rm = TRUE),
  sd = ~sd(.x, na.rm = TRUE),
  median = ~median(.x, na.rm = TRUE),
  n.avail = ~sum(!is.na(.x))
)
