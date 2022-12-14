% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/co_net.R
\name{c_net_cal}
\alias{c_net_cal}
\title{Calculate spearman correlation for one or two t(otutab)}
\usage{
c_net_cal(totu, totu2 = NULL, filename = "occor", threads = 4)
}
\arguments{
\item{totu}{t(otutab)}

\item{totu2}{t(otutab) or NULL}

\item{filename}{the prefix of saved files}

\item{threads}{parallel mode threads}

\item{parallel}{open parallel mode?(default:F)}
}
\value{
a list with 2 elements:
\item{r}{spearman correlation}
\item{p.value}{p.adjust.method = 'NULL}
}
\description{
Calculate spearman correlation for one or two t(otutab)
}
\examples{
data("otutab")
t(otutab) -> totu
c_net_cal(totu) -> corr
metadata[,12:19] -> env
c_net_cal(totu, env) -> corr2
}
