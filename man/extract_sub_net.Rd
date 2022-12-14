% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/co_net.R
\name{extract_sub_net}
\alias{extract_sub_net}
\title{Extract sub-nwtwork from the whole network}
\usage{
extract_sub_net(a_net, otutab, threads = 4, save_net = F)
}
\arguments{
\item{a_net}{the whole network}

\item{otutab}{otutab, these columns will be extract}

\item{threads}{threads, default:4}

\item{save_net}{should save these sub_nets? F or a filename}
}
\value{
a dataframe contains all sub_net parameters
}
\description{
Extract sub-nwtwork from the whole network
}
\examples{
data(otutab)
c_net_cal(t(otutab))->corr
corr$p.value<-p.adjust.table(corr$p.value)
c_net_build(corr)->a_net
extract_sub_net(a_net,otutab,save_net = "testnet")
}
