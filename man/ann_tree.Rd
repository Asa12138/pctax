% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pc_phylogenetic.R
\name{ann_tree}
\alias{ann_tree}
\alias{easy_tree}
\title{Annotate a tree}
\usage{
ann_tree(f_tax, otutab, level = 7)

easy_tree(tree)
}
\arguments{
\item{f_tax}{taxonomy dataframe}

\item{otutab}{otutab, rowname==colname(taxonomy)}

\item{level}{1~7}

\item{tree}{result from ann_tree}
}
\value{
a treedata

a ggplot
}
\description{
Annotate a tree

Easy way to plot a phylogenetic tree
}
\examples{
data(otutab)
ann_tree(taxonomy,otutab)->tree
easy_tree(tree)
}
