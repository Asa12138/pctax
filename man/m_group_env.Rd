% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional.R
\name{m_group_env}
\alias{m_group_env}
\alias{plot.mant_g}
\title{Multi-table test with env}
\usage{
m_group_env(g_otutab, env)

\method{plot}{mant_g}(x, ...)
}
\arguments{
\item{g_otutab}{multi-otutabs with first column is group}

\item{env}{environmental factors}

\item{x}{mant_g object}

\item{...}{add}
}
\value{
a mant_g object

a ggplot
}
\description{
Multi-table test with env

Plot mant_g object
}
\examples{
if (requireNamespace("linkET")) {
  data(otutab, package = "pcutils")
  cbind(group = rep(c("a", "b", "c"), c(200, 100, 185)), otutab) -> g_otutab
  metadata[, 3:8, drop = FALSE] -> env
  m_group_env(g_otutab, env) -> mant_g
  plot(mant_g)
}
}
