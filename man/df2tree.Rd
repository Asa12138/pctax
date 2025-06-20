% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogenetic.R
\name{df2tree}
\alias{df2tree}
\title{From a dataframe to construct a phylo}
\usage{
df2tree(data, edge_df = FALSE, ignore_pattern = NULL)
}
\arguments{
\item{data}{dataframe}

\item{edge_df}{if the data is edge_df ?}

\item{ignore_pattern}{An optional regular expression pattern to match tip or node labels for dropping.}
}
\value{
phylo object
}
\description{
NOTE: this function will do \code{before_tree} first.
}
\examples{
data(otutab, package = "pcutils")
df2tree(taxonomy) -> tax_tree
print(tax_tree)
# check all nodes matched!
if (requireNamespace("picante")) {
  picante::match.phylo.comm(tax_tree, t(otutab)) -> nn
  nrow(nn$comm) == nrow(t(otutab))
}
}
