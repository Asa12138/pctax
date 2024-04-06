#' Create a pc_otu class object
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param metadata a metadata data.frame, samples are rows
#' @param taxonomy a taxomomy data.frame, look out the rowname of taxonomy and otutab should matched!
#' @param ... add
#'
#' @export
#' @return pc_otu
#' @examples
#' data(otutab, package = "pcutils")
#' pc_tax1 <- pc_otu(otutab, metadata)
#' print(pc_tax1)
pc_otu <- function(otutab = data.frame(), metadata = data.frame(), taxonomy = NULL, ...) {
  if (!is.data.frame(otutab)) stop("otutab must be a df!")
  if (!is.data.frame(metadata)) stop("metadata must be a df!")
  if (!is.null(taxonomy)) {
    if (!is.data.frame(taxonomy)) stop("taxonomy must be a df!")
    if (length(taxonomy) != 7) stop("taxonomy should be a df with 'k,p,c,o,f,g,s'!")
    pre_tax_table(taxonomy) -> taxonomy
    taxonomy[match(rownames(otutab), rownames(taxonomy)), ] %>% na.omit() -> taxonomy
  }
  ids <- intersect(rownames(metadata), colnames(otutab))
  if (length(ids) < 1) stop("no samples left, check the rownames of metadata or colnames of otutable")
  metadata <- metadata[ids, , drop = FALSE]
  otutab <- otutab[, ids, drop = FALSE]
  otutab <- otutab[rowSums(otutab) > 0, , drop = FALSE]
  pc <- list(
    tbls = list(otutab = otutab),
    metas = list(metadata = metadata),
    otus = list(taxonomy = taxonomy), ...
  )
  class(pc) <- append("pc_otu", class(pc))
  return(pc)
}

#' Judge pc_otu is valid or not
#'
#' @param pc a pc_otu object
#' @return logical
#' @export
pc_valid <- function(pc) {
  for (i in pc$tbls) if (!is.data.frame(i) & !is.null(i)) stop("tbls must be df!")
  for (i in pc$metas) if (!is.data.frame(i) & !is.null(i)) stop("metas must be df!")
  for (i in pc$otus) if (!is.data.frame(i) & !is.null(i)) stop("otus must be df!")
  if (!all(rownames(pc$metas$metadata) == colnames(pc$tbls$otutab))) stop("Wrong samples! check the rownames of metadata or colnames of otutable")
  return(TRUE)
}

#' Print
#'
#' @param x pc_otu
#' @param ... add
#'
#' @return No value
#' @method print pc_otu
#' @exportS3Method
print.pc_otu <- function(x, ...) {
  pc <- x
  sprintf("There are %d otus and %d samples!", nrow(pc$tbls$otutab), ncol(pc$tbls$otutab)) %>% print()
  for (i in names(pc)) {
    pcutils::dabiao(i, print = TRUE)
    if (i %in% c("tbls", "metas", "otus")) {
      for (j in names(pc[[i]])) {
        pcutils::dabiao(j, n = 40, print = TRUE)
        if ("data.frame" %in% class(pc[[i]][[j]])) {
          print(head(pc[[i]][[j]]))
        } else {
          print(pc[[i]][[j]])
        }
      }
    } else {
      if ("data.frame" %in% class(pc[[i]])) {
        print(head(pc[[i]]))
      } else {
        print(pc[[i]])
      }
    }
  }
}

#' Summary pc_otu
#' @param object pc_otu
#' @param ... add
#' @return No value
#' @method summary pc_otu
#' @exportS3Method
#' @examples
#' data("pc_tax1")
#' summary(pc_tax1)
#'
summary.pc_otu <- function(object, ...) {
  pc <- object
  pc_valid(pc)
  sprintf("There are %d otus and %d samples!", nrow(pc$tbls$otutab), ncol(pc$tbls$otutab)) %>% print()
  pcutils::dabiao("tables, include some data filter and tranformat", print = TRUE)
  print(names(pc$tbls))
  pcutils::dabiao("metadatas, include some statistics", print = TRUE)
  print(names(pc$metas))
  pcutils::dabiao("otus annotation, include some statistics", print = TRUE)
  print(names(pc$otus))
  pcutils::dabiao("some other indexs:", print = TRUE)
  print(names(pc)[-1:-3])
}

#' Add taxonomy for a pc_otu object
#' @export
#' @param pc a pc_otu object
#' @param taxonomy a taxomomy data.frame, look out the rownames of taxonomy and otutab should matched!
#' @return pc_otu
#' @examples
#' data(otutab, package = "pcutils")
#' pc_tax1 <- pc_otu(otutab, metadata)
#' pc_tax1 <- add_tax(pc_tax1, taxonomy)
add_tax <- function(pc, taxonomy) {
  if (!"pc_otu" %in% class(pc)) stop(pc, "should be a pc_otu")
  if (!is.data.frame(taxonomy)) stop("taxonomy must be a df!")
  if (length(taxonomy) != 7) stop("taxonomy should be a df with 'k,p,c,o,f,g,s'!")

  pre_tax_table(taxonomy) -> taxonomy
  taxonomy[match(rownames(pc$tbls$otutab), rownames(taxonomy)), ] %>% na.omit() -> taxonomy
  pc$otus$taxonomy <- taxonomy
  return(pc)
}
