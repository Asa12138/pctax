# phylogenic tree==============

#' Complete a taxonomy table
#'
#' @param tax_table taxonomy table
#' @param tax_levels a vector whose length longer than `ncol(taxdf)`, use to be prefix. Default: c("k", "p", "c", "o", "f", "g","s", "st")
#' @param na_tax grepl some words and turn to `na_repalce`, default: "Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig"
#' @param ignore.case ignore.case for `na_tax`
#' @param na_repalce defalut: Unknown
#'
#' @return a good taxonomy table
#' @export
#' @references \code{MicrobiotaProcess}
#' @examples
#' taxmat <- matrix(sample("onelevel", 7 * 2, replace = TRUE), nrow = 2, ncol = 7) %>% as.data.frame()
#' colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' pre_tax_table(taxmat)
pre_tax_table <- function(tax_table, tax_levels = c("k", "p", "c", "o", "f", "g", "s", "st"),
                          na_tax = "Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig", ignore.case = TRUE,
                          na_repalce = "Unknown") {
  if (length(tax_levels) < ncol(tax_table)) stop("need tax_levels length at least: ", ncol(tax_table))
  if (any(duplicated(tax_levels))) stop("tax_levels can not be duplicated.")

  rown <- rownames(tax_table)
  # 找到部分匹配的tax变成NA，如果一列全是NA则应该去掉
  tax_table <- as.data.frame(tax_table, check.names = FALSE)
  if (is.na(na_tax) | is.null(na_tax)) na_tax <- ""
  if (nchar(na_tax) > 0) tax_table[pcutils::grepl.data.frame(na_tax, tax_table, ignore.case = ignore.case)] <- NA
  indx <- vapply(tax_table, function(i) all(is.na(i)), FUN.VALUE = logical(1)) %>% setNames(NULL)
  tax_table <- tax_table[, !indx, drop = FALSE]

  # 添加tax_levels
  if (!(grepl(paste0("^", tax_levels[1], "__"), tax_table[1, 1]))) {
    tmprownames <- rownames(tax_table)
    tmpcolnames <- colnames(tax_table)
    tax_table <- t(apply(tax_table, 1, as.character))
    tax_table[is.na(tax_table)] <- ""
    tax_table <- data.frame(t(apply(tax_table, 1, \(i)paste(tax_levels[seq_len(length(i))], i, sep = "__"))),
      stringsAsFactors = FALSE
    )
    rownames(tax_table) <- tmprownames
    colnames(tax_table) <- tmpcolnames
  }

  # 把NA变成Unknown
  tax_table1 <- tax_table
  tmprownames <- rownames(tax_table1)
  # about 2chars more
  indexmark <- apply(tax_table1, 2, \(x) {
    nchar(x, keepNA = TRUE)
  }) <= 3
  tax_table1[indexmark] <- NA

  if (any(is.na(tax_table1[, 1]))) {
    prefix <- paste0(tax_levels[1], "__")
    tax_table1[is.na(tax_table1[, 1]), 1] <- paste0(prefix, na_repalce)
  }

  indextmp <- apply(is.na(tax_table1), 1, which)
  if (length(indextmp) == 0) {
    tax_table1 <- data.frame(tax_table1, check.names = FALSE)
    return(tax_table1)
  }

  tax_table1 <- apply(tax_table1, 1, zoo::na.locf)
  tax_table1 <- lapply(seq_len(ncol(tax_table1)), function(i) tax_table1[, i])
  # 如果本身是NA，就应该变为形如g__un_f__onelevel
  tax_table1 <- data.frame(t(mapply(function(x, y) {
    y <- as.vector(y)
    x[y] <- paste(tax_levels[y], x[y], sep = "__un_")
    x
  }, tax_table1, indextmp)), stringsAsFactors = FALSE)
  rownames(tax_table1) <- tmprownames

  # 查找有没有一些tax的parent不同导致的重复，有的话就要变成sametax_A,sametax_B
  for (i in seq_len(7)) {
    tax_table2 <- tax_table1
    if (ncol(tax_table2) == 1) {
      return(tax_table2)
    }
    tax_table2 <- tax_table2 %>% rownames_to_column()
    for (i in ncol(tax_table2):3) {
      tmp <- split(tax_table2, tax_table2[, i])
      for (j in seq_len(length(tmp))) {
        flag <- length(unique(as.vector(tmp[[j]][, i - 1])))
        if (flag > 1) {
          tmp[[j]][, i] <- paste(tmp[[j]][, i], tmp[[j]][, i - 1], sep = "_")
        }
      }
      tax_table2 <- do.call("rbind", c(tmp, make.row.names = FALSE))
    }
    tax_table1 <- tax_table2 %>% column_to_rownames(var = "rowname")
  }

  tax_table <- tax_table1[rown, ]
  return(tax_table)
}


#' Calculate the lowest common ancestor (LCA) of a set of taxa
#'
#' @param df a data frame with taxonomic information, with columns representing taxonomic levels
#'
#' @return character
#' @export
#'
#' @examples
#' df <- data.frame(
#'   A = c("a", "a", "a", "a"),
#'   B = c("x", "x", "y", "y"),
#'   C = c("1", "1", "2", "3"),
#'   stringsAsFactors = FALSE
#' )
#' tax_lca(df)
tax_lca <- function(df) {
  # 检查输入是否为数据框
  if (!is.data.frame(df)) {
    stop("Input must be a data frame.")
  }

  # 获取列数
  num_cols <- ncol(df)

  # 从最后一列开始向前遍历
  for (i in num_cols:1) {
    # 获取当前列的唯一值
    unique_values <- unique(df[[i]])

    # 如果唯一值的数量为1，则返回该值
    if (length(unique_values) == 1) {
      return(unique_values)
    }
  }

  # 如果所有列都没有相同的值，返回NA或其他标志
  return("_root")
}

#' Before df2tree check
#'
#' @param f_tax table
#'
#' @return table
#' @export
#'
#' @examples
#' wrong_taxdf <- data.frame(
#'   kingdom = c(rep(c("A", "B"), each = 4), "C", NA),
#'   "phylum" = c("A", "a", "b", "c", "c", "c", "d", NA, NA, "e")
#' )
#' before_tree(wrong_taxdf)
#'
before_tree <- function(f_tax) {
  f_tax <- as.data.frame(f_tax)

  du_name <- c()
  exist_name <- unique(f_tax[, 1]) %>% na.omit()
  for (i in 2:ncol(f_tax)) {
    tmp <- unique(f_tax[, i]) %>% na.omit()
    tmp <- intersect(tmp, exist_name)
    if (length(tmp) > 0) {
      du_name <- c(du_name, tmp)
      ind <- f_tax[, i] %in% tmp
      # f_tax[ind,i]=paste0(f_tax[ind,i],strrep(" ", i - 1))
      f_tax[ind, i] <- paste0(colnames(f_tax)[i], "_", f_tax[ind, i])
    }
    exist_name <- c(exist_name, unique(f_tax[, i]))
  }

  if (length(du_name) > 0) message("Some name exist in different column:\n", paste(du_name, collapse = ", "))

  na_parent <- c()
  for (i in seq_len(ncol(f_tax) - 1)) {
    tmp <- dplyr::distinct(f_tax[, c(i, i + 1)])
    idx <- is.na(tmp[, 1, drop = T]) & !is.na(tmp[, 2, drop = T])
    du <- tmp[idx, 2]
    if (length(du) > 0) {
      na_parent <- c(na_parent, paste0(colnames(f_tax)[i + 1], ":", du))
      for (idx1 in which(idx)) {
        f_tax[, i] <- apply(f_tax, 1, \(x)ifelse(identical(unlist(tmp[idx1, ]), x[c(i, i + 1)]), paste(x[i + 1], "parent", sep = "_"), x[i]))
      }
    }
  }

  if (length(na_parent) > 0) message("Some name have NA parent:\n", paste(na_parent, collapse = ", "))

  diff_parent_name <- c()

  for (i in seq_len(ncol(f_tax) - 1)) {
    tmp <- dplyr::distinct(f_tax[, c(i, i + 1)])
    du <- tmp[duplicated(tmp[, 2]), 2] %>% na.omit()

    if (length(du) > 0) {
      diff_parent_name <- c(diff_parent_name, du)
      f_tax[, i + 1] <- apply(f_tax, 1, \(x)ifelse(x[i + 1] %in% du, paste(x[i + 1], x[i], sep = "_"), x[i + 1]))
    }
    is.na(tmp[, 1, drop = T]) & !is.na(tmp[, 2, drop = T])
  }
  if (length(diff_parent_name) > 0) message("Some name have different parents:\n", paste(diff_parent_name, collapse = ", "))

  f_tax
}


#' From a dataframe to construct a phylo
#' @description
#' NOTE: this function will do `before_tree` first.
#'
#' @param data dataframe
#' @param edge_df if the data is edge_df ?
#'
#' @return phylo object
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' df2tree(taxonomy) -> tax_tree
#' print(tax_tree)
#' # check all nodes matched!
#' if (requireNamespace("picante")) {
#'   picante::match.phylo.comm(tax_tree, t(otutab)) -> nn
#'   nrow(nn$comm) == nrow(t(otutab))
#' }
df2tree <- function(data, edge_df = FALSE) {
  if (!edge_df) {
    data <- before_tree(data)

    data <- data.frame(Root = rep("r__root", nrow(data)), data)
    datalist <- list()
    clnm <- colnames(data)
    for (i in seq_len(ncol(data) - 1)) {
      tmpdat <- data[, c(i, i + 1)]
      colnames(tmpdat) <- c("parent", "child")
      # tmpdat %<>% dplyr::mutate(nodeClass = clnm[i + 1], nodeDepth = i) %>%
      #     dplyr::distinct()
      datalist[[i]] <- tmpdat
    }
    datalist <- do.call("rbind", datalist)
  } else {
    datalist <- data
    colnames(datalist) <- c("parent", "child")
  }


  datalist <- datalist[!duplicated(datalist), ]
  isTip <- !as.vector(datalist$child) %in% as.vector(datalist$parent)
  index <- rep(NA, length(isTip))
  index[isTip] <- seq(1, sum(isTip))
  index[!isTip] <- seq(sum(isTip) + 2, length(isTip) + 1)
  mapping <- data.frame(
    node = index, labelnames = as.vector(datalist$child),
    isTip
  )
  indxx <- match(mapping$labelnames, datalist$child)
  # mapping$nodeClass <- datalist[indxx, "nodeClass"]
  # mapping$nodeDepth <- datalist[indxx, "nodeDepth"]
  parentnode <- mapping[match(
    as.vector(datalist$parent),
    as.vector(mapping$labelnames)
  ), ]$node
  childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)), ]$node
  edges <- cbind(parentnode, childnode)
  colnames(edges) <- NULL
  edges[is.na(edges)] <- sum(isTip) + 1
  root <- data.frame(
    node = sum(isTip) + 1, labelnames = "r__root",
    isTip = FALSE
    # nodeClass = "Root", nodeDepth = 0
  )
  mapping <- rbind(root, mapping)
  mapping <- mapping[order(mapping$node), ]

  node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
  tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
  # mapping <- mapping[, colnames(mapping) %in% c(
  #     "node",
  #     "nodeClass", "nodeDepth"
  # )]
  taxphylo <- structure(
    list(
      edge = edges, node.label = node.label, edge.length = rep(1, nrow(edges)),
      tip.label = tip.label, Nnode = length(node.label)
    ),
    class = "phylo"
  )
  return(taxphylo)
}

# 会将带空格的label改变，淘汰

#' From a dataframe to construct a phylo (save nwk)
#' @description
#' NOTE: this function will transfer all space to `_`
#'
#' @param taxa dataframe
#' @return phylo object
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' df2tree(taxonomy) -> tax_tree
#' print(tax_tree)
df2tree1 <- function(taxa) {
  # taxa%>%mutate_all(.funs = \(x)gsub(" ","",x))->taxa
  makeNewick1 <- \(taxa, naSub = "_") {
    if (!is.null(naSub)) {
      taxa[is.na(taxa)] <- naSub
    }
    if (ncol(taxa) == 0) {
      return("")
    }
    bases <- unique(taxa[, 1])
    innerTree <- sapply(bases, function(ii) makeNewick1(taxa[taxa[, 1] == ii, -1, drop = FALSE]))
    out <- sprintf("(%s)", paste(sprintf("%s%s", innerTree, bases), collapse = ","))
    return(out)
  }
  taxa <- pcutils::gsub.data.frame(" ", "_", taxa)
  paste0(makeNewick1(taxa), ";") -> tree
  tree %>% ape::read.tree(text = .) -> nwk
  nwk$edge.length <- rep(1, nrow(nwk$edge))
  return(nwk)
}

#' Annotate a tree
#'
#' @param f_tax taxonomy dataframe
#' @param otutab otutab, rowname==rowname(taxonomy)
#' @param level 1~7
#'
#' @return a treedata
#' @export
#'
#' @examples
#' if (interactive()) {
#'   data(otutab, package = "pcutils")
#'   ann_tree(taxonomy, otutab) -> tree
#'   # run yourself
#'   easy_tree(tree, add_abundance = FALSE) -> p
#'   p
#' }
ann_tree <- function(f_tax, otutab = NULL, level = ncol(f_tax)) {
  lib_ps("ggtree", "vctrs", library = FALSE)
  label.x <- label.y <- NULL
  # le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  le <- c("root", colnames(f_tax))
  df2tree(f_tax[, 1:level]) %>% ggtree::fortify() -> tree

  tree$level <- le[ceiling(tree$branch) + 1]
  # tree$level<-factor(tree$level,levels = le)
  # dplyr操作后会把tbl_tree class删除
  dplyr::left_join(tree, tree[, c("node", "label")], by = c("parent" = "node")) %>%
    dplyr::rename(label = label.x, parent_label = label.y) -> tree1

  if (!is.null(otutab)) {
    # if(any(rownames(f_tax)!=rownames(otutab)))stop("rowname not match")
    otutab <- otutab[rownames(f_tax), , drop = FALSE]
  } else {
    otutab <- data.frame(row.names = rownames(f_tax), n = rep(1, nrow(f_tax)))
  }
  otutab %>% rowSums(.) / ncol(otutab) -> num
  res <- data.frame(label = "", abundance = sum(num))
  for (i in 1:level) {
    aggregate(num, by = list(f_tax[, i]), sum) -> tmp
    colnames(tmp) <- colnames(res)
    res <- rbind(res, tmp)
  }
  dplyr::left_join(tree1, res, by = "label") -> tree2

  ann_tax <- data.frame()
  for (i in level:1) {
    f_tax[, 1:i, drop = FALSE] %>% dplyr::distinct() -> tmpdf
    tmpdf[, ncol(tmpdf)] -> rownames(tmpdf)
    ann_tax <- vctrs::vec_c(ann_tax, tmpdf)
  }
  dplyr::left_join(tree2, ann_tax %>% tibble::rownames_to_column("label"), by = "label") -> tree3
  if (!"tbl_tree" %in% class(tree3)) class(tree3) <- c("tbl_tree", class(tree3))
  return(tree3)
}

#' Easy way to plot a phylogenetic tree
#'
#' @param tree result from `ann_tree`
#' @param highlight highlight which level, one of `tree$level`
#' @param colorfill "color" or "fill"
#' @param pal color pal
#' @param add_abundance logical
#' @param add_tiplab logical
#' @param basic_params parameters parse to \code{\link[ggtree]{ggtree}}
#' @param fontsize tip label fontsize
#' @param name_prefix keep the prefix like "k__" or "p__" in the label? Default: FALSE
#' @param color_name color name
#' @param topN topN to show
#'
#' @importFrom ggplot2 geom_tile
#' @return a ggplot
#' @export
#'
#' @rdname ann_tree
easy_tree <- function(tree, highlight = "Phylum", colorfill = "color", topN = NULL, pal = NULL, name_prefix = FALSE,
                      basic_params = NULL, add_abundance = TRUE, color_name = "abundance", add_tiplab = TRUE, fontsize = NULL) {
  label <- level <- node <- in_label <- group <- abundance <- phy_label <- isTip <- NULL
  lib_ps("ggtree", "ggtreeExtra", library = FALSE)
  requireNamespace("ggplot2")
  colorfill <- match.arg(colorfill, c("color", "fill"))

  # 可以剪切部分树
  if (!name_prefix) {
    tree %>% dplyr::mutate(`in_label` = sub("^.?__", "", label)) -> tree2
  } else {
    tree %>% dplyr::mutate(`in_label` = label) -> tree2
  }

  dplyr::filter(dplyr::as_tibble(tree2), level == highlight) %>% dplyr::select(node, `in_label`) -> phy
  colnames(phy)[2] <- "phy_label"

  if (!is.null(topN)) {
    if (is.numeric(topN) && topN > 0) {
      tree2 %>%
        dplyr::filter(isTip) %>%
        dplyr::count(.data[[highlight]]) -> tmp_count
      dplyr::top_n(tmp_count, n = topN, n) %>% dplyr::pull(highlight) -> topN
      if (!name_prefix) {
        topN <- sub("^.?__", "", topN)
      }
    }
    phy %>% dplyr::mutate(`phy_label` = ifelse(phy_label %in% topN, phy_label, "Others")) -> phy
  }
  n_phy <- unique(phy$phy_label) %>% length()

  if (is.null(fontsize)) {
    fontsize <- round(600 / sum(tree2$isTip), 2)
    if (fontsize > 5) fontsize <- 5
    if (fontsize < 0.2) fontsize <- 0.2
  } else if (!is.numeric(fontsize)) stop("fontsize should be numeric")

  default_pal <- setNames(pcutils::get_cols(n_phy), unique(phy$phy_label))

  if (!is.null(pal)) {
    if (is.null(names(pal))) {
      names(pal) <- rep(unique(phy$phy_label), length = length(pal))
    }
  }

  la <- setdiff(unique(phy$phy_label), names(pal))
  if (length(la) > 0) message("Some labels not in names(pal): ", paste(la, collapse = ", "), "\nSo set the default color.")

  pal <- pcutils::update_param(default_pal, pal)

  if (colorfill == "fill") {
    p <- do.call(ggtree::ggtree, pcutils::update_param(list(tr = tree2, layout = "radial", size = 0.5), basic_params)) +
      ggtree::geom_highlight(data = phy, aes(node = node, fill = `phy_label`), alpha = 0.6, to.bottom = TRUE) +
      scale_fill_manual(name = highlight, values = pal, na.value = NA)
  }

  if (colorfill == "color") {
    tree2 <- ggtree::groupClade(tree2, setNames(phy$node, phy$`phy_label`))
    tree2$group <- as.character(tree2$group)
    tree2$group[tree2$group == 0] <- "_root"

    legends <- levels(as.factor(tree2$group) %>% pcutils::change_fac_lev("_root", last = TRUE))
    tree2$group <- factor(tree2$group, levels = legends)

    basic_color <- "black"
    basic_linetype <- 1
    if ("color" %in% names(basic_params)) basic_color <- ifelse(pcutils::is.ggplot.color(basic_params$color), basic_params$color, "black")
    if ("col" %in% names(basic_params)) basic_color <- ifelse(pcutils::is.ggplot.color(basic_params$col), basic_params$col, "black")

    if (is.null(names(pal))) {
      pal <- setNames(c(pal, basic_color), legends)
    } else {
      pal <- c(pal, "_root" = basic_color)
    }

    if ("linetype" %in% names(basic_params)) basic_linetype <- ifelse(is.numeric(basic_params$linetype), basic_params$linetype, 1)

    p <- do.call(ggtree::ggtree, pcutils::update_param(list(tr = tree2, mapping = aes(color = group), layout = "radial", size = 0.5), basic_params)) +
      scale_color_manual(name = highlight, values = pal, na.translate = FALSE, label = c(legends[-length(legends)], "")) +
      # geom_point()+
      guides(colour = guide_legend(
        order = 1,
        override.aes = list(
          linewidth = 2,
          linetype = setNames(c(rep(basic_linetype, n_phy), NA), legends)
        )
      ))
  }
  if (add_abundance & ("abundance" %in% colnames(tree2))) {
    geom_tile <- ggplot2::geom_tile
    p <- p + ggnewscale::new_scale_fill() +
      ggtreeExtra::geom_fruit(
        geom = geom_tile,
        mapping = aes(y = label, fill = abundance),
        stat = "identity", offset = 0.08, pwidth = 0.1,
      ) + scale_fill_gradient(name = color_name, low = "#FFFFFF", high = "#B2182B")
  }
  if (add_tiplab) {
    p <- p + ggtree::geom_tiplab(aes(label = `in_label`), color = "black", size = fontsize, offset = ifelse(add_abundance, 1.1, 0.2), show.legend = FALSE)
  }
  p
}


# Get strip for a circle tree
get_strip <- \(name, tree, flat_n = 5){
  label <- node <- x <- NULL
  lib_ps("tidytree", library = FALSE)
  stopifnot(inherits(tree, "tbl_tree"))
  mx <- max(tree$x)
  id <- dplyr::filter(tree, label == name) %>% dplyr::pull(node)
  tidytree::offspring(tree, id) %>%
    dplyr::filter(x == mx) %>%
    dplyr::arrange(angle) -> tmp
  offset.text <- ifelse(nrow(tmp) > flat_n, 0.5, 0)
  hjust <- ifelse(nrow(tmp) > flat_n, 0.5, 0)
  angle <- ifelse(nrow(tmp) > flat_n, mean(tmp$angle) - 90, mean(tmp$angle))
  list(head(tmp, 1) %>% dplyr::pull(label), tail(tmp, 1) %>% dplyr::pull(label), offset.text, hjust, angle)
}

#' add strips for a tree plot
#'
#' @param trp tree plot from `ggtree`
#' @param some_tax some tax you want to add strip
#' @param flat_n flat the text when taxa number more than `flat_n`.
#' @param strip_params parameters parse to \code{\link[ggtree]{geom_strip}}
#'
#' @return tree plot
#' @export
#'
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' # run yourself
#' if (interactive()) {
#'   ann_tree(taxonomy, otutab) -> tree
#'   easy_tree(tree) -> p
#'   some_tax <- table(taxonomy$Phylum) %>%
#'     sort(decreasing = TRUE) %>%
#'     head(5) %>%
#'     names()
#'   add_strip(p, some_tax)
#' }
#' }
add_strip <- function(trp, some_tax, flat_n = 5, strip_params = NULL) {
  trp1 <- trp
  tree <- trp$data

  # add strip
  # some_tax=table(read_taxonomy$Phylum)%>%sort(decreasing = TRUE)%>%head(10)%>%names()

  cols <- pcutils::get_cols(length(some_tax), c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"))

  for (i in seq_len(length(some_tax))) {
    add_f <- get_strip(some_tax[i], tree, flat_n = flat_n)
    trp1 <- trp1 +
      do.call(ggtree::geom_strip, pcutils::update_param(list(
        taxa1 = add_f[[1]], taxa2 = add_f[[2]],
        barsize = 2, color = cols[i], label = gsub("^.?__", "", some_tax[i]), fontsize = 3,
        offset = 1, offset.text = add_f[[3]], hjust = add_f[[4]], angle = add_f[[5]]
      ), strip_params))
  }
  trp1
}


#' Plot two trees in one plot
#'
#' @param tree1 phylo object
#' @param tree2 phylo object
#' @param edge_df dataframe with edge information, containing "from" and "to" columns
#' @param tree2_x x position of tree2
#' @param filter_link filter the link between tree1 and tree2
#' @param tree1_param parameters for \code{\link[ggtree]{geom_tree}}
#' @param tree2_param parameters for \code{\link[ggtree]{geom_tree}}
#' @param line_param parameters for \code{\link[ggplot2]{geom_line}}
#' @param tree1_tip tree tip label
#' @param tip1_param parameters for \code{\link[ggtree]{geom_tiplab}}
#' @param tree2_tip tree tip label
#' @param tip2_param parameters for \code{\link[ggtree]{geom_tiplab}}
#' @param tree1_highlight tree1 highlight data.frame
#' @param highlight1_param parameters for \code{\link[ggtree]{geom_hilight}}
#' @param highlight1_scale scale_fill_ for highlight1
#' @param tree2_highlight tree2 highlight data.frame
#' @param highlight2_param parameters for \code{\link[ggtree]{geom_hilight}}
#' @param highlight2_scale scale_fill_ for highlight2
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' df2tree(taxonomy[1:50, ]) -> tax_tree
#' df2tree(taxonomy[51:100, ]) -> tax_tree2
#' link <- data.frame(from = sample(tax_tree$tip.label, 20), to = sample(tax_tree2$tip.label, 20))
#' plot_two_tree(tax_tree, tax_tree2, link)
plot_two_tree <- function(tree1, tree2, edge_df = NULL, tree2_x = 10, filter_link = FALSE,
                          tree1_param = list(), tree2_param = list(),
                          line_param = list(),
                          tree1_tip = FALSE, tip1_param = list(), tree2_tip = FALSE, tip2_param = list(),
                          tree1_highlight = NULL, highlight1_param = list(), highlight1_scale = NULL,
                          tree2_highlight = NULL, highlight2_param = list(), highlight2_scale = ggplot2::scale_fill_hue(na.value = NA)) {
  if (!is.null(edge_df)) {
    if (!all(c("from", "to") %in% colnames(edge_df))) {
      stop("edge_df must have columns 'from' and 'to'")
    }
    edge_df <- mutate_if(edge_df, is.factor, as.character)
    if (!"id" %in% colnames(edge_df)) {
      edge_df$id <- 1:nrow(edge_df)
    } else if (anyDuplicated(edge_df$id)) {
      stop("edge_df must have unique 'id' column")
    }
    if (!"width" %in% colnames(edge_df)) edge_df$width <- 1
    if (!"e_type" %in% colnames(edge_df)) edge_df$e_type <- "correlation"

    if (filter_link) {
      tree1 <- ape::drop.tip(tree1, setdiff(tree1$tip.label, edge_df$from))
      if (is.null(tree1)) stop("No tips left in tree1")
      tree2 <- ape::drop.tip(tree2, setdiff(tree2$tip.label, edge_df$to))
      if (is.null(tree2)) stop("No tips left in tree2")
    }
  }

  # tree
  tr1 <- do.call(ggtree::ggtree, update_param(list(tr = tree1), tree1_param))$data
  tr2 <- do.call(ggtree::ggtree, update_param(list(tr = tree2), tree1_param))$data

  tr2$x <- max(tr2$x) - tr2$x + max(tr1$x) + tree2_x
  tr1$y <- pcutils::mmscale(tr1$y, min(tr2$y), max(tr2$y))

  if (!is.null(edge_df)) {
    trdf1 <- filter(tr1, !is.na(label), isTip) %>% data.frame(row.names = .$label, .)
    trdf2 <- filter(tr2, !is.na(label), isTip) %>% data.frame(row.names = .$label, .)

    lapply(seq_len(nrow(edge_df)), \(x){
      i <- edge_df[x, ]
      data.frame(
        rbind(trdf1[i$from, ], trdf2[i$to, ]),
        group = i$id,
        width = i$width,
        e_type = i$e_type
      )
    }) %>% do.call(rbind, .) -> tmp_edge
  }

  p <- do.call(ggtree::ggtree, update_param(list(tr = tr1), tree1_param))

  if (!is.null(tree1_highlight)) {
    if (!all(c("id", "group") %in% colnames(tree1_highlight))) stop("tree1_highlight must have columns 'node' and 'group'")
    p <- p + do.call(ggtree::geom_hilight, update_param(
      list(
        data = tree1_highlight %>% select(id, group),
        mapping = aes(node = id, fill = group)
      ),
      highlight1_param
    ))
    if (!is.null(highlight1_scale)) p <- p + highlight1_scale
  }

  p <- p + do.call(ggtree::geom_tree, update_param(list(data = tr2), tree2_param))

  if (!is.null(tree2_highlight)) {
    if (!all(c("id", "group") %in% colnames(tree2_highlight))) stop("tree2_highlight must have columns 'id' and 'group'")
    if (!is.null(tree1_highlight)) p <- p + ggnewscale::new_scale_fill()
    p <- p + do.call(ggtree::geom_hilight, update_param(
      list(
        data = left_join(tr2, tree2_highlight, by = c("node" = "id")),
        mapping = aes(node = node, fill = group)
      ),
      highlight2_param
    ))
    if (!is.null(highlight2_scale)) p <- p + highlight2_scale
  }

  if (!is.null(edge_df)) {
    p <- p + do.call(
      ggplot2::geom_line,
      update_param(
        list(mapping = aes(x, y, group = group, color = e_type, linewidth = width), data = tmp_edge),
        line_param
      )
    )
  }

  if (tree1_tip) {
    p <- p + do.call(ggtree::geom_tiplab, update_param(list(geom = "text"), tip1_param))
  }
  if (tree2_tip) {
    p <- p + do.call(ggtree::geom_tiplab, update_param(list(data = tr2, hjust = 1), tip2_param))
  }
  p + scale_linewidth_continuous(range = c(0.5, 2))
}
