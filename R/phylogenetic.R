# phylogenic tree==============
#
taxaclass <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Rank8", "Rank9", "Rank10")

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
    lib_ps("dplyr", "tibble", library = FALSE)
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
    tax_table1=tax_table
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
        tax_table2=tax_table1
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

#' Before df2tree check
#'
#' @param f_tax table
#'
#' @return table
#' @export
#'
#' @examples
#' wrong_taxdf <- data.frame(
#'     kingdom = c(rep(c("A", "B"), each = 4), "C", NA),
#'     "phylum" = c("A", "a", "b", "c", "c", "c", "d", NA, NA, "e")
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
#' @param data dataframe
#'
#' @return phylo object
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' df2tree(taxonomy) -> tax_tree
#' print(tax_tree)
#' # check all nodes matched!
#' picante::match.phylo.comm(tax_tree, t(otutab)) -> nn
#' nrow(nn$comm) == nrow(t(otutab))
df2tree <- function(data) {
    data <- before_tree(data)

    data <- data.frame(Root = rep("r__root", nrow(data)), data)
    datalist <- list()
    clnm <- colnames(data)
    for (i in seq_len(ncol(data) - 1)) {
        tmpdat <- data[, c(i, i + 1)]
        colnames(tmpdat) <- c("parent", "child")
        tmpdat %<>% dplyr::mutate(nodeClass = clnm[i + 1], nodeDepth = i) %>%
            dplyr::distinct()
        datalist[[i]] <- tmpdat
    }
    datalist <- do.call("rbind", datalist)
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
    mapping$nodeClass <- datalist[indxx, "nodeClass"]
    mapping$nodeDepth <- datalist[indxx, "nodeDepth"]
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
        isTip = FALSE, nodeClass = "Root", nodeDepth = 0
    )
    mapping <- rbind(root, mapping)
    mapping <- mapping[order(mapping$node), ]
    node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
    tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
    mapping <- mapping[, colnames(mapping) %in% c(
        "node", "nodeClass",
        "nodeDepth"
    )]
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
#' data(otutab, package = "pcutils")
#' ann_tree(taxonomy, otutab) -> tree
#' # run yourself
#' easy_tree(tree, add_abundance = FALSE) -> p
#' p
#'
ann_tree <- function(f_tax, otutab = NULL, level = ncol(f_tax)) {
    lib_ps("ggtree", "vctrs", library = FALSE)

    # le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    le <- c("root", colnames(f_tax))
    df2tree(f_tax[, 1:level]) %>% ggtree::fortify() -> tree

    tree$level <- le[ceiling(tree$branch) + 1]
    # tree$level<-factor(tree$level,levels = le)
    dplyr::left_join(tree, tree[, c("node", "label")], by = c("parent" = "node")) %>% dplyr::rename(label = label.x, parent_label = label.y) -> tree1

    if (!is.null(otutab)) {
        # if(any(rownames(f_tax)!=rownames(otutab)))stop("rowname not match")
        otutab <- otutab[rownames(f_tax), , drop = FALSE]
        otutab %>% rowSums(.) / ncol(otutab) -> num

        res <- data.frame(label = "", abundance = sum(num))
        for (i in 1:level) {
            aggregate(num, by = list(f_tax[, i]), sum) -> tmp
            colnames(tmp) <- colnames(res)
            res <- rbind(res, tmp)
        }
        dplyr::left_join(tree1, res) -> tree2
    } else {
        tree2 <- tree1
    }

    ann_tax <- data.frame()
    for (i in level:1) {
        f_tax[, 1:i, drop = FALSE] %>% dplyr::distinct() -> tmpdf
        tmpdf[, ncol(tmpdf)] -> rownames(tmpdf)
        ann_tax <- vctrs::vec_c(ann_tax, tmpdf)
    }
    dplyr::left_join(tree2, ann_tax %>% tibble::rownames_to_column("label"), by = "label") -> tree3

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
#'
#' @importFrom ggplot2 geom_tile
#' @return a ggplot
#' @export
#'
#' @rdname ann_tree
easy_tree <- function(tree, highlight = "Phylum", colorfill = "color", pal = NULL, basic_params = NULL, add_abundance = TRUE, add_tiplab = TRUE, fontsize = NULL) {
    lib_ps("ggtree", "ggtreeExtra", library = FALSE)
    requireNamespace("ggplot2")
    colorfill <- match.arg(colorfill, c("color", "fill"))

    # 可以剪切部分树
    tree %>% dplyr::mutate(`internal_label` = sub("^.?__", "", label)) -> tree2
    dplyr::filter(tree2, level == highlight) %>% dplyr::select(node, `internal_label`) -> phy

    if (is.null(fontsize)) {
        fontsize <- round(600 / sum(tree2$isTip), 2)
        if (fontsize > 5) fontsize <- 5
        if (fontsize < 0.2) fontsize <- 0.2
    } else if (!is.numeric(fontsize)) stop("fontsize should be numeric")

    if (!is.null(pal)) {
        if (length(pal) < nrow(phy)) stop("need ", nrow(phy), " colors, just give ", length(pal))
    } else {
        pal <- pcutils::get_cols(nrow(phy))
    }

    if (colorfill == "fill") {
        p <- do.call(ggtree::ggtree, pcutils::update_param(list(tr = tree2, layout = "radial", size = 0.5), basic_params)) +
            ggtree::geom_highlight(data = phy, aes(node = node, fill = `internal_label`), alpha = 0.6, to.bottom = TRUE) +
            scale_fill_manual(name = highlight, values = pal, na.value = NA)
    }

    if (colorfill == "color") {
        tree2 <- ggtree::groupClade(tree2, setNames(phy$node, phy$`internal_label`))
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
                    linetype = setNames(c(rep(basic_linetype, nrow(phy)), NA), legends)
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
            ) + scale_fill_gradient(name = "abundance", low = "#FFFFFF", high = "#B2182B")
    }
    if (add_tiplab) {
        p <- p + ggtree::geom_tiplab(aes(label = `internal_label`), color = "black", size = fontsize, offset = ifelse(add_abundance, 1.1, 0.2), show.legend = FALSE)
    }
    p
}


# Get strip for a circle tree
get_strip <- \(name, tree, flat_n = 5){
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
#' ann_tree(taxonomy, otutab) -> tree
#' # run yourself
#' if (FALSE) {
#'     easy_tree(tree) -> p
#'     some_tax <- table(taxonomy$Phylum) %>%
#'         sort(decreasing = TRUE) %>%
#'         head(5) %>%
#'         names()
#'     add_strip(p, some_tax)
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


#' Plot a sankey
#'
#' @param tree result from \code{\link{ann_tree}}
#' @param top_N each level has top_N
#' @param notshow some words you don't want to show
#' @param intermediate logical, show the intermediate rank
#' @param width width
#' @param height height
#' @param ... look for parameters in \code{\link[sankeyD3]{sankeyNetwork}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' ann_tree(taxonomy[, c(1, 5, 6, 7)], otutab) -> tree
#' sangji_plot(tree)
#' }
sangji_plot <- function(tree, top_N = 5, notshow = c(), intermediate = FALSE, width = 3000, height = 500, ...) {
    lib_ps("sankeyD3", library = FALSE)
    lib_ps("tidytree", library = FALSE)

    # 桑基图
    le <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    # select show part
    if (length(notshow) > 0) {
        tree1 <- tree[!grepl(paste0(notshow, collapse = "|"), tree$label), ]
    } else {
        tree1 <- tree
    }

    tree1 %>%
        dplyr::group_by(level) %>%
        dplyr::top_n(top_N, abundance) %>%
        dplyr::ungroup() -> sangji_dat

    if (intermediate) {
        all_node <- sapply(sangji_dat$node, \(i)tidytree::ancestor(tree, i) %>% dplyr::pull(node)) %>% unlist()
        all_node <- c(sangji_dat$node, all_node) %>% unique()
        tree %>% dplyr::filter(node %in% all_node) -> sangji_dat
    }

    # show levels
    sangji_dat %>%
        dplyr::filter(!label %in% c("", "r__root")) %>%
        dplyr::select(label, level, abundance) %>%
        as.data.frame() -> nodes

    # tree structure
    sangji_dat %>%
        dplyr::filter(parent_label != "r__root") %>%
        dplyr::select(label, parent_label, abundance) -> links
    links <- as.data.frame(links)

    # find the nearest parent
    others <- links$parent_label[!links$parent_label %in% nodes$label]
    others <- (others[!duplicated(others)])
    tree %>% dplyr::filter(label %in% others) -> o_nodes
    others <- o_nodes$label
    o_nodes <- o_nodes$node
    for (i in seq_along(o_nodes)) {
        tidytree::ancestor(tree, o_nodes[i]) %>% dplyr::pull(label) -> tmp
        links[links$parent_label == others[i], "parent_label"] <- rev(tmp)[rev(tmp) %in% nodes$label][1]
    }

    le[le %in% nodes$level] -> mytax

    taxRank_to_depth <- setNames(seq_along(mytax) - 1, mytax)
    nodes$depth <- taxRank_to_depth[nodes$level %>% as.character()]


    links$IDsource <- match(links$parent_label, nodes$label) - 1
    links$IDtarget <- match(links$label, nodes$label) - 1
    # na.omit(links)->links

    do.call(sankeyD3::sankeyNetwork, pcutils::update_param(list(
        Links = links, Nodes = nodes,
        Source = "IDsource", Target = "IDtarget", Value = "abundance",
        NodeID = "label", NodeGroup = "label", NodePosX = "depth", NodeValue = "abundance",
        iterations = 1000, xAxisDomain = mytax, align = "none",
        fontSize = 12, linkGradient = TRUE,
        nodeWidth = 15, nodeCornerRadius = 5, highlightChildLinks = TRUE,
        orderByPath = TRUE, scaleNodeBreadthsByString = TRUE,
        numberFormat = "pavian", dragY = TRUE, nodeShadow = TRUE,
        doubleclickTogglesChildren = TRUE, width = width, height = height
    ), list(...)))
}

#' Plot a sunburst
#'
#' @param tree result from \code{\link{ann_tree}}
#'
#' @return sunburst
#' @export
#' @seealso [sangji_plot()]
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' ann_tree(taxonomy[, c(1, 5, 6, 7)], otutab) -> tree
#' sunburst(tree)
#' }
sunburst <- function(tree) {
    lib_ps("plotly", library = FALSE)
    links <- data.frame("source" = tree$parent_label, "target" = tree$label, "weight" = tree$abundance)

    fig <- plotly::plot_ly(
        labels = links$target, parents = links$source,
        values = links$weight, type = "sunburst"
    )
    fig
}
