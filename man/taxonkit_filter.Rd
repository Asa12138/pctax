% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rtaxonkit.R
\name{taxonkit_filter}
\alias{taxonkit_filter}
\title{Filter TaxIDs based on Taxonomic Ranks}
\usage{
taxonkit_filter(
  file_path,
  black_list = NULL,
  discard_noranks = FALSE,
  discard_root = FALSE,
  equal_to = NULL,
  higher_than = NULL,
  lower_than = NULL,
  rank_file = NULL,
  root_taxid = NULL,
  save_predictable_norank = FALSE,
  taxid_field = NULL,
  text = FALSE,
  data_dir = NULL
)
}
\arguments{
\item{file_path}{The path to the input file containing TaxIDs. Or file text (text=TRUE)}

\item{black_list}{A character vector specifying the ranks to discard.}

\item{discard_noranks}{Logical value indicating whether to discard all ranks without order (default is FALSE).}

\item{discard_root}{Logical value indicating whether to discard the root taxid (default is FALSE).}

\item{equal_to}{A character vector specifying the ranks for which TaxIDs should be equal to.}

\item{higher_than}{The rank above which the TaxIDs should be (exclusive).}

\item{lower_than}{The rank below which the TaxIDs should be (exclusive).}

\item{rank_file}{The path to a user-defined ordered taxonomic ranks file.}

\item{root_taxid}{The root taxid (default is 1).}

\item{save_predictable_norank}{Logical value indicating whether to save some special ranks without order when using lower_than (default is FALSE).}

\item{taxid_field}{The field index of the taxid in the input file (default is 1).}

\item{text}{logical}

\item{data_dir}{directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")}
}
\value{
A character vector containing the output of the "taxonkit filter" command.
}
\description{
This function uses the "taxonkit filter" command to filter TaxIDs based on taxonomic ranks.
}
\examples{
\dontrun{
taxids2 <- system.file("extdata/taxids2.txt", package = "pctax")
taxonkit_filter(taxids2, lower_than = "genus")
}
}
\seealso{
Other Rtaxonkit: 
\code{\link{check_taxonkit}()},
\code{\link{download_taxonkit_dataset}()},
\code{\link{install_taxonkit}()},
\code{\link{name_or_id2df}()},
\code{\link{taxonkit_lca}()},
\code{\link{taxonkit_lineage}()},
\code{\link{taxonkit_list}()},
\code{\link{taxonkit_name2taxid}()},
\code{\link{taxonkit_reformat}()}
}
\concept{Rtaxonkit}
