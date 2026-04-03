#' @keywords internal
"_PACKAGE"


## usethis namespace: start
#' @importFrom utils data
#' @importFrom utils head tail
#' @importFrom stats aggregate median var sd setNames runif relevel time na.omit kmeans anova t.test
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom pcutils lib_ps hebing trans dabiao generate_labels update_param get_cols t2 change_fac_lev
#' @import dplyr
#' @import ggplot2
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "id"))

pctax_theme <- {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      plot.margin = grid::unit(rep(0.5, 4), "lines"),
      strip.background = ggplot2::element_rect(fill = NA)
    )
}
