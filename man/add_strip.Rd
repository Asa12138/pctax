% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogenetic.R
\name{add_strip}
\alias{add_strip}
\title{add strips for a tree plot}
\usage{
add_strip(trp, some_tax, flat_n = 5, strip_params = NULL)
}
\arguments{
\item{trp}{tree plot from \code{ggtree}}

\item{some_tax}{some tax you want to add strip}

\item{flat_n}{flat the text when taxa number more than \code{flat_n}.}

\item{strip_params}{parameters parse to \code{\link[ggtree]{geom_strip}}}
}
\value{
tree plot
}
\description{
add strips for a tree plot
}
\examples{
\donttest{
data(otutab, package = "pcutils")
# run yourself
if (interactive()) {
  ann_tree(taxonomy, otutab) -> tree
  easy_tree(tree) -> p
  some_tax <- table(taxonomy$Phylum) \%>\%
    sort(decreasing = TRUE) \%>\%
    head(5) \%>\%
    names()
  add_strip(p, some_tax)
}
}
}
