% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a_diversity.R
\name{z_diversity_decay}
\alias{z_diversity_decay}
\alias{plot.zeta_decay}
\title{Calculate Zeta Diversity with Distance}
\usage{
z_diversity_decay(otutab, xy_df, group_df = NULL, zetadiv_params = list())

\method{plot}{zeta_decay}(x, ribbon = TRUE, ...)
}
\arguments{
\item{otutab}{A matrix or data frame containing OTU (Operational Taxonomic Unit) counts.}

\item{xy_df}{Site coordinates.}

\item{group_df}{A data frame containing group information.}

\item{zetadiv_params}{Additional parameters to be passed to the Zeta.ddecay function from the zetadiv package.}

\item{x}{Zeta diversity results obtained from z_diversity_decay function.}

\item{ribbon}{Logical, whether to add a ribbon to the plot for standard deviation.}

\item{...}{Additional arguments to be passed to ggplot2 functions.}
}
\value{
zeta_decay

A ggplot object.
}
\description{
This function calculates Zeta diversity for each group in the provided otutab.
}
\examples{
if (requireNamespace("zetadiv")) {
  data(otutab, package = "pcutils")
  zeta_decay_result <- z_diversity_decay(otutab, metadata[, c("lat", "long")],
    metadata["Group"],
    zetadiv_params = list(sam = 10)
  )
  plot(zeta_decay_result)
}
}
