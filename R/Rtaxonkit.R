#' Install taxonkit
#'
#' @param taxonkit_tar_gz your download taxonkit_tar_gz file from https://github.com/shenwei356/taxonkit/releases/
#' @param make_sure make sure to do this
#' @family Rtaxonkit
#' @return No value
#' @export
install_taxonkit <- function(make_sure = FALSE, taxonkit_tar_gz = NULL) {
  # Detect the operating system
  os <- tolower(Sys.info()[["sysname"]])
  machine <- ifelse(grepl("arm", tolower(Sys.info()[["machine"]])), "arm",
    ifelse(grepl("amd", tolower(Sys.info()[["machine"]])), "amd", "others")
  )

  # Set the installation URL and command based on the operating system
  if (os == "linux") {
    url <- "https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_linux_amd64.tar.gz"
    if (machine == "arm") url <- "https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_linux_arm64.tar.gz"
    install_cmd <- "tar -xf taxonkit_linux_*.tar.gz"
  } else if (os == "darwin") {
    url <- "https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_darwin_amd64.tar.gz"
    if (machine == "arm") url <- "https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_darwin_arm64.tar.gz"
    install_cmd <- "tar -xf taxonkit_darwin_*.tar.gz"
  } else if (os == "windows") {
    url <- "https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_windows_amd64.exe.tar.gz"
  } else {
    stop("Unsupported operating system!")
  }
  # Set the file name and installation directory
  filename <- basename(url)

  install_dir <- tools::R_user_dir("pctax")
  install_dir <- normalizePath(install_dir) %>% suppressWarnings()

  if (!make_sure) {
    message("please set `make_sure=TRUE` to install taxonkit at ", install_dir)
    return(invisible())
  }
  # Create the installation directory if it does not exist
  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)

  # Set the full path to the installation directory
  install_path <- file.path(install_dir, filename)

  if (is.null(taxonkit_tar_gz)) {
    # Download the file
    pcutils::download2(url, install_path)
  } else {
    if (file.exists(taxonkit_tar_gz) & grepl("taxonkit_.*\\.tar\\.gz", taxonkit_tar_gz)) {
      file.copy(taxonkit_tar_gz, install_path)
    } else {
      stop("Wrong file: ", taxonkit_tar_gz)
    }
  }

  # Extract the downloaded file
  ori_dir <- getwd()
  on.exit(setwd(ori_dir))

  setwd(install_dir) # Change to the installation directory
  if (os == "windows") {
    utils::untar("taxonkit_windows_amd64.exe.tar.gz")
  } else {
    system(install_cmd)
  }

  # Remove the downloaded file
  file.remove(install_path)

  message("taxonkit has been successfully installed!\n")
}

#' Download taxonkit dataset
#'
#' @param taxdump_tar_gz your download taxdump_tar_gz file from https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
#' @param make_sure make sure to do this
#' @family Rtaxonkit
#' @export
#' @return No value
download_taxonkit_dataset <- function(make_sure = FALSE, taxdump_tar_gz = NULL) {
  url <- "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  home_dir <- Sys.getenv("HOME")
  dest_dir <- file.path(home_dir, ".taxonkit")

  taxdump <- file.path(tools::R_user_dir("pctax"), "taxdump.tar.gz")
  taxdump <- normalizePath(taxdump) %>% suppressWarnings()

  if (!make_sure) {
    message("please set `make_sure=TRUE` to download taxonkit dataset at ", taxdump)
    return(invisible())
  }
  # Download the taxdump.tar.gz file
  if (is.null(taxdump_tar_gz)) {
    ori_time <- getOption("timeout")
    on.exit(options(timeout = ori_time))

    options(timeout = 300)
    tryCatch(expr = {
      utils::download.file(url, destfile = taxdump, mode = "wb")
    }, error = function(e) {
      stop("Try download yourself from https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz ")
    })
    # pcutils::download2(url,taxdump,mode = "wb")
  } else {
    if (file.exists(taxdump_tar_gz) & grepl("taxdump.tar.gz", taxdump_tar_gz)) {
      taxdump <- taxdump_tar_gz
    } else {
      stop("Wrong file: ", taxdump_tar_gz)
    }
  }

  if (file.exists(file.path(dest_dir, "names.dmp"))) {
    message("Taxonkit dataset already exists at ", dest_dir, "\nReplace it?")
    if (!utils::askYesNo()) {
      return(invisible())
    }
  }

  # Uncompress the tar.gz file
  utils::untar(taxdump, exdir = dest_dir)

  # Copy the required files to the destination directory
  # files_to_copy <- c("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp")
  # for (file in files_to_copy) {
  #   file_path <- file.path("taxdump", file)
  #   dest_path <- file.path(dest_dir, file)
  #   file.copy(file_path, dest_path)
  # }

  # Clean up temporary files
  # file.remove("taxdump.tar.gz")
  # unlink(taxdump, recursive = TRUE)

  message("Taxonkit files downloaded and copied successfully.\n")
}

#' Check taxonkit
#'
#' @param print print
#' @family Rtaxonkit
#' @export
#' @return taxonkit path
check_taxonkit <- function(print = TRUE) {
  taxonkit <- file.path(tools::R_user_dir("pctax"), "taxonkit")
  taxonkit <- normalizePath(taxonkit) %>% suppressWarnings()
  if (!file.exists(taxonkit)) stop("Taxonkit not found, please try `install_taxonkit()`")
  flag <- system(paste(shQuote(taxonkit), "-h"), ignore.stdout = !print, ignore.stderr = TRUE)
  if (flag != 0) stop("Taxonkit not found, please try `install_taxonkit()`")
  if (print) pcutils::dabiao("Taxonkit is available if there is help message above")

  home_dir <- Sys.getenv("HOME")
  dest_dir <- file.path(home_dir, ".taxonkit")
  if (dir.exists(dest_dir) & all(c("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp") %in% list.files(dest_dir))) {
    if (print) pcutils::dabiao("Taxonkit dataset is available!")
  } else {
    stop("Taxonkit dataset (", dest_dir, ") not found, please try `download_taxonkit_dataset()`")
  }
  if (!print) {
    return(taxonkit)
  }
}

#' Taxonkit list
#'
#' This function uses Taxonkit to perform the "list" operation, which retrieves
#' information about taxa based on their TaxIDs.
#'
#' @param ids A character vector of TaxIDs to retrieve information for.
#' @param indent The indentation string to use for pretty-printing the output. Default is "  ".
#' @param json Logical value indicating whether to output the result in JSON format. Default is FALSE.
#' @param show_name Logical value indicating whether to show the scientific names of taxa. Default is FALSE.
#' @param show_rank Logical value indicating whether to show the ranks of taxa. Default is FALSE.
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return The output of the Taxonkit list operation.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' taxonkit_list(ids = c(9605), indent = "-", show_name = TRUE, show_rank = TRUE)
#' }
taxonkit_list <- function(ids, indent = "  ", json = FALSE, show_name = FALSE, show_rank = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)
  # Convert the ids to a comma-separated string
  ids_str <- paste(ids, collapse = ",")

  # Build the taxonkit command
  cmd <- paste(shQuote(taxonkit), "list")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (json) cmd <- paste(cmd, "-J")
  if (show_name) cmd <- paste(cmd, "-n")
  if (show_rank) cmd <- paste(cmd, "-r")
  cmd <- paste(cmd, "-i", ids_str)
  if (indent != "  ") cmd <- paste(cmd, " -I '", indent, "'", sep = "")

  # Execute the taxonkit command
  res <- system(cmd, intern = TRUE)
  res
}


#' Retrieve Taxonomic Lineage using taxonkit
#'
#' @param file_path The path to the input file with taxonomic IDs. Or file text (text=TRUE)
#' @param delimiter The field delimiter in the lineage (default ";").
#' @param no_lineage Logical, indicating whether to exclude lineage information (default: FALSE).
#' @param show_lineage_ranks Logical, indicating whether to append ranks of all levels in the lineage (default: FALSE).
#' @param show_lineage_taxids Logical, indicating whether to append lineage consisting of taxids (default: FALSE).
#' @param show_name Logical, indicating whether to append scientific name (default: FALSE).
#' @param show_rank Logical, indicating whether to append rank of taxids (default: FALSE).
#' @param show_status_code Logical, indicating whether to show status code before lineage (default: FALSE).
#' @param taxid_field The field index of taxid. Input data should be tab-separated (default: 1).
#' @param text logical,
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return A character vector containing the taxonomic lineage information.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' taxonkit_lineage("9606\n63221", show_name = TRUE, show_rank = TRUE, text = TRUE)
#' }
taxonkit_lineage <- function(file_path, delimiter = ";", no_lineage = FALSE, show_lineage_ranks = FALSE,
                             show_lineage_taxids = FALSE, show_name = FALSE, show_rank = FALSE,
                             show_status_code = FALSE, taxid_field = 1, text = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)

  # Prepare taxonkit command
  taxonkit_cmd <- paste(shQuote(taxonkit), "lineage")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (text) {
    tmp_path <- file.path(tempdir(), "Rtaxonkit.tmp")
    cat(file_path, sep = "\n", file = tmp_path)
    file_path <- tmp_path
  }

  # Add options based on user inputs
  if (no_lineage) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-L", sep = " ")
  }
  if (show_lineage_ranks) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-R", sep = " ")
  }
  if (show_lineage_taxids) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-t", sep = " ")
  }
  if (show_name) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-n", sep = " ")
  }
  if (show_rank) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-r", sep = " ")
  }
  if (show_status_code) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-c", sep = " ")
  }
  if (delimiter != ";") {
    taxonkit_cmd <- paste(taxonkit_cmd, "-d", delimiter)
  }
  if (taxid_field != 1) {
    taxonkit_cmd <- paste(taxonkit_cmd, "-i", taxid_field)
  }

  # Execute taxonkit command on the input file
  cmd <- paste(taxonkit_cmd, file_path, sep = " ")
  lineage <- system(cmd, intern = TRUE)

  # Return the lineage
  lineage
}

#' Reformat Taxonomic Lineage using taxonkit
#'
#' @param file_path The path to the input file with taxonomic lineages. Or file text (text=TRUE)
#' @param delimiter The field delimiter in the input lineage (default ";").
#' @param add_prefix Logical, indicating whether to add prefixes for all ranks (default: FALSE).
#' @param prefix_kingdom The prefix for kingdom, used along with --add-prefix (default: "K__").
#' @param prefix_phylum The prefix for phylum, used along with --add-prefix (default: "p__").
#' @param prefix_class The prefix for class, used along with --add-prefix (default: "c__").
#' @param prefix_order The prefix for order, used along with --add-prefix (default: "o__").
#' @param prefix_family The prefix for family, used along with --add-prefix (default: "f__").
#' @param prefix_genus The prefix for genus, used along with --add-prefix (default: "g__").
#' @param prefix_species The prefix for species, used along with --add-prefix (default: "s__").
#' @param prefix_subspecies The prefix for subspecies, used along with --add-prefix (default: "t__").
#' @param prefix_strain The prefix for strain, used along with --add-prefix (default: "T__").
#' @param fill_miss_rank Logical, indicating whether to fill missing rank with lineage information of the next higher rank (default: FALSE).
#' @param format_string The output format string with placeholders for each rank.
#' @param miss_rank_repl_prefix The prefix for estimated taxon level for missing rank (default: "unclassified ").
#' @param miss_rank_repl The replacement string for missing rank.
#' @param miss_taxid_repl The replacement string for missing taxid.
#' @param output_ambiguous_result Logical, indicating whether to output one of the ambiguous result (default: FALSE).
#' @param lineage_field The field index of lineage. Input data should be tab-separated (default: 2).
#' @param taxid_field The field index of taxid. Input data should be tab-separated. It overrides -i/--lineage-field.
#' @param pseudo_strain Logical, indicating whether to use the node with lowest rank as strain name (default: FALSE).
#' @param trim Logical, indicating whether to not fill missing rank lower than current rank (default: FALSE).
#' @param text logical
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return A character vector containing the reformatted taxonomic lineages.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' # Use taxid
#' taxids2 <- system.file("extdata/taxids2.txt", package = "pctax")
#' reformatted_lineages <- taxonkit_reformat(taxids2,
#'   add_prefix = TRUE, taxid_field = 1, fill_miss_rank = TRUE
#' )
#' reformatted_lineages
#' taxonomy <- strsplit2(reformatted_lineages, "\t")
#' taxonomy <- strsplit2(taxonomy$V2, ";")
#'
#' # Use lineage result
#' taxonkit_lineage("9606\n63221", show_name = TRUE, show_rank = TRUE, text = TRUE) %>%
#'   taxonkit_reformat(text = TRUE)
#' }
taxonkit_reformat <- function(file_path,
                              delimiter = NULL,
                              add_prefix = FALSE,
                              prefix_kingdom = "K__",
                              prefix_phylum = "p__",
                              prefix_class = "c__",
                              prefix_order = "o__",
                              prefix_family = "f__",
                              prefix_genus = "g__",
                              prefix_species = "s__",
                              prefix_subspecies = "t__",
                              prefix_strain = "T__",
                              fill_miss_rank = FALSE,
                              format_string = "",
                              miss_rank_repl_prefix = "unclassified ",
                              miss_rank_repl = "",
                              miss_taxid_repl = "",
                              output_ambiguous_result = FALSE,
                              lineage_field = 2,
                              taxid_field = NULL,
                              pseudo_strain = FALSE,
                              trim = FALSE, text = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)
  # Prepare taxonkit command
  cmd <- paste(shQuote(taxonkit), "reformat")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (text) {
    tmp_path <- file.path(tempdir(), "Rtaxonkit.tmp")
    cat(file_path, sep = "\n", file = tmp_path)
    file_path <- tmp_path
  }

  if (add_prefix) {
    cmd <- paste0(
      cmd, " --add-prefix",
      " --prefix-K ", prefix_kingdom,
      " --prefix-p ", prefix_phylum,
      " --prefix-c ", prefix_class,
      " --prefix-o ", prefix_order,
      " --prefix-f ", prefix_family,
      " --prefix-g ", prefix_genus,
      " --prefix-s ", prefix_species,
      " --prefix-t ", prefix_subspecies,
      " --prefix-T ", prefix_strain
    )
  }
  if (fill_miss_rank) {
    cmd <- paste0(cmd, " --fill-miss-rank")
  }
  if (format_string != "") {
    cmd <- paste0(cmd, " --format '", format_string, "'")
  }
  if (miss_rank_repl != "") {
    cmd <- paste0(cmd, " --miss-rank-repl '", miss_rank_repl, "'")
  }
  if (miss_rank_repl_prefix != "unclassified ") {
    cmd <- paste0(cmd, " --miss-rank-repl-prefix '", miss_rank_repl_prefix, "'")
  }
  if (miss_taxid_repl != "") {
    cmd <- paste0(cmd, " --miss-taxid-repl '", miss_taxid_repl, "'")
  }
  if (output_ambiguous_result) {
    cmd <- paste0(cmd, " --output-ambiguous-result")
  }
  cmd <- paste0(cmd, " --lineage-field ", lineage_field)
  if (!is.null(taxid_field)) {
    cmd <- paste0(cmd, " --taxid-field ", taxid_field)
  }
  if (pseudo_strain) {
    cmd <- paste0(cmd, " --pseudo-strain")
  }
  if (trim) {
    cmd <- paste0(cmd, " --trim")
  }
  if (!is.null(delimiter)) {
    cmd <- paste0(cmd, " -d '", delimiter, "'")
  }

  # Execute the command
  reformatted_lineages <- system(paste(cmd, file_path), intern = TRUE)

  # Return the reformatted lineages
  reformatted_lineages
}

#' Convert Taxonomic Names to TaxIDs
#'
#' This function uses the "taxonkit taxonkit_name2taxid" command to convert taxonomic names to corresponding taxonomic IDs (TaxIDs).
#'
#' @param file_path The path to the input file containing taxonomic names. Or file text (text=TRUE)
#' @param name_field The field index of the taxonomic name in the input file (default is 1).
#' @param sci_name Logical value indicating whether to search only for scientific names (default is FALSE).
#' @param text Logical
#' @param show_rank Logical value indicating whether to show the taxonomic rank in the output (default is FALSE).
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return A character vector containing the output of the "taxonkit_name2taxid" command.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' names <- system.file("extdata/name.txt", package = "pctax")
#' taxonkit_name2taxid(names, name_field = 1, sci_name = FALSE, show_rank = FALSE)
#' "Homo sapiens" %>% taxonkit_name2taxid(text = TRUE)
#' }
taxonkit_name2taxid <- function(file_path, name_field = NULL, sci_name = FALSE, show_rank = FALSE, text = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)

  # Prepare taxonkit command
  cmd <- paste(shQuote(taxonkit), "name2taxid")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (text) {
    tmp_path <- file.path(tempdir(), "Rtaxonkit.tmp")
    cat(file_path, sep = "\n", file = tmp_path)
    file_path <- tmp_path
  }

  if (!is.null(name_field)) {
    cmd <- paste0(cmd, " --name-field ", name_field)
  }
  if (sci_name) {
    cmd <- paste0(cmd, " --sci-name")
  }
  if (show_rank) {
    cmd <- paste0(cmd, " --show-rank")
  }

  # Execute the command
  output <- system(paste(cmd, file_path), intern = TRUE)

  # Return the output
  output
}

#' Filter TaxIDs based on Taxonomic Ranks
#'
#' This function uses the "taxonkit filter" command to filter TaxIDs based on taxonomic ranks.
#'
#' @param file_path The path to the input file containing TaxIDs. Or file text (text=TRUE)
#' @param black_list A character vector specifying the ranks to discard.
#' @param discard_noranks Logical value indicating whether to discard all ranks without order (default is FALSE).
#' @param discard_root Logical value indicating whether to discard the root taxid (default is FALSE).
#' @param equal_to A character vector specifying the ranks for which TaxIDs should be equal to.
#' @param higher_than The rank above which the TaxIDs should be (exclusive).
#' @param lower_than The rank below which the TaxIDs should be (exclusive).
#' @param rank_file The path to a user-defined ordered taxonomic ranks file.
#' @param root_taxid The root taxid (default is 1).
#' @param save_predictable_norank Logical value indicating whether to save some special ranks without order when using lower_than (default is FALSE).
#' @param text logical
#' @param taxid_field The field index of the taxid in the input file (default is 1).
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return A character vector containing the output of the "taxonkit filter" command.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' taxids2 <- system.file("extdata/taxids2.txt", package = "pctax")
#' taxonkit_filter(taxids2, lower_than = "genus")
#' }
taxonkit_filter <- function(file_path, black_list = NULL, discard_noranks = FALSE, discard_root = FALSE,
                            equal_to = NULL, higher_than = NULL, lower_than = NULL, rank_file = NULL,
                            root_taxid = NULL, save_predictable_norank = FALSE, taxid_field = NULL, text = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)
  # Prepare taxonkit command
  cmd <- paste(shQuote(taxonkit), "filter")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (text) {
    tmp_path <- file.path(tempdir(), "Rtaxonkit.tmp")
    cat(file_path, sep = "\n", file = tmp_path)
    file_path <- tmp_path
  }

  if (!is.null(black_list)) {
    cmd <- paste0(cmd, " -B '", paste(black_list, collapse = ","), "'")
  }
  if (discard_noranks) {
    cmd <- paste0(cmd, " -N")
  }
  if (discard_root) {
    cmd <- paste0(cmd, " -R")
  }
  if (!is.null(equal_to)) {
    cmd <- paste0(cmd, " -E '", paste(equal_to, collapse = ","), "'")
  }
  if (!is.null(higher_than)) {
    cmd <- paste0(cmd, " -H '", higher_than, "'")
  }
  if (!is.null(lower_than)) {
    cmd <- paste0(cmd, " -L '", lower_than, "'")
  }
  if (!is.null(rank_file)) {
    cmd <- paste0(cmd, " -r '", rank_file, "'")
  }
  if (!is.null(root_taxid)) {
    cmd <- paste0(cmd, " --root-taxid ", root_taxid)
  }
  if (save_predictable_norank) {
    cmd <- paste0(cmd, " -n")
  }
  if (!is.null(taxid_field)) {
    cmd <- paste0(cmd, " -i ", taxid_field)
  }

  # Execute the command
  output <- system(paste(cmd, file_path), intern = TRUE)

  # Return the output
  output
}

#' Compute Lowest Common Ancestor (LCA) of TaxIDs
#'
#' This function uses the "taxonkit lca" command to compute the Lowest Common Ancestor (LCA) of TaxIDs.
#'
#' @param file_path The path to the input file containing TaxIDs. Or file text (text=TRUE)
#' @param buffer_size The size of the line buffer (supported units: K, M, G).
#' @param separator The separator for TaxIDs.
#' @param skip_deleted Whether to skip deleted TaxIDs and compute with the remaining ones.
#' @param skip_unfound Whether to skip unfound TaxIDs and compute with the remaining ones.
#' @param taxids_field The field index of TaxIDs. Input data should be tab-separated (default 1).
#' @param text logical
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#'
#' @return A character vector containing the computed LCAs.
#' @export
#' @family Rtaxonkit
#' @examples
#' \dontrun{
#' taxonkit_lca("239934, 239935, 349741", text = TRUE, separator = ", ")
#' }
taxonkit_lca <- function(file_path, buffer_size = "1M", separator = " ",
                         skip_deleted = FALSE, skip_unfound = FALSE, taxids_field = NULL, text = FALSE, data_dir = NULL) {
  taxonkit <- check_taxonkit(print = FALSE)

  # Prepare taxonkit command
  cmd <- paste(shQuote(taxonkit), "lca")
  if (!is.null(data_dir)) cmd <- paste(cmd, " --data-dir ", data_dir, sep = "")

  if (text) {
    tmp_path <- file.path(tempdir(), "Rtaxonkit.tmp")
    cat(file_path, sep = "\n", file = tmp_path)
    file_path <- tmp_path
  }

  command <- paste(
    cmd,
    "--buffer-size", buffer_size,
    "--separator", separator,
    if (skip_deleted) "--skip-deleted" else "",
    if (skip_unfound) "--skip-unfound" else "",
    if (!is.null(taxids_field)) paste("--taxids-field", taxids_field) else "",
    file_path
  )
  # Execute the command
  output <- system(command, intern = TRUE)

  # Return the output
  output
}

#' Transfer taxon name or taxid to the lineage dataframe
#'
#' @param name_or_id name or taxid
#' @param mode "id" or "name"
#' @param data_dir directory containing nodes.dmp and names.dmp (default "/Users/asa/.taxonkit")
#' @param add_prefix add_prefix
#' @param fill_miss_rank fill_miss_rank
#' @family Rtaxonkit
#' @return dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' name_or_id2df(c("Homo sapiens", "Akkermansia muciniphila ATCC BAA-835"))
#' }
name_or_id2df <- function(name_or_id, mode = "name", add_prefix = TRUE, fill_miss_rank = TRUE, data_dir = NULL) {
  if (mode == "name") {
    df <- name_or_id %>%
      taxonkit_name2taxid(text = TRUE, data_dir = data_dir) %>%
      utils::read.table(text = ., sep = "\t", col.names = c("name", "taxid"), comment.char = "")
  } else if (mode == "id") df <- data.frame(taxid = name_or_id)
  reformatted_lineages <- taxonkit_reformat(df$taxid,
    add_prefix = add_prefix, fill_miss_rank = fill_miss_rank, text = TRUE,
    taxid_field = 1, data_dir = data_dir
  )
  taxonomy <- pcutils::strsplit2(reformatted_lineages, "\t")
  taxonomy <- pcutils::strsplit2(taxonomy$V2, ";", colnames = c(
    "Kingdom", "Phylum", "Class", "Order", "Family",
    "Genus", "Species"
  ))
  cbind(df, taxonomy)
}
