## some suggested pkgs: start
# AnnotationDbi,
# org.Hs.eg.db,
# linkET,
# pairwiseAdonis,
# sankeyD3,
# ggchicklet,
# Biobase,
# GEOquery
## some suggested pkgs: end

#' Permanova between a otutab and a variable (added `two`)
#'
#' @param two two by two adonis test
#' @inheritParams permanova
permanova2 <- function(otutab, envs, norm = TRUE, each = TRUE, method = "adonis",
                       dist = "bray", two = FALSE, nperm = 999, ...) {
  all <- c("adonis", "anosim", "mrpp", "mantel")
  if (!method %in% all) stop(paste0("method should be one of ", paste0(all, collapse = ", ")))
  stopifnot(is.data.frame(envs))

  match_res <- match_df(otutab, envs)
  otutab <- match_res$otutab
  env <- match_res$metadata

  data.frame(t(otutab)) -> dat
  if (norm) otu.t <- vegan::decostand(dat, "hellinger") else otu.t <- dat
  if (each) {
    # adnois不只是检验分组变量，连续的数值变量也可以检验。
    soil <- NULL
    if (method == "adonis") {
      for (i in 1:ncol(env)) {
        dat.div <- vegan::adonis2(otu.t ~ (env[, i]), permutations = nperm, method = dist)
        soil <- rbind(soil, c(colnames(env)[i], dat.div$R2[1], dat.div$`Pr(>F)`[1]))
        if (two) {
          if ((is.factor(env[, i]) | inherits(env[, i], "Date") | is.character(env[, i]))) {
            env[, i] %>% as.factor() -> group
            lib_ps("pairwiseAdonis", library = FALSE)
            dat.pairwise.adonis <- pairwiseAdonis::pairwise.adonis(
              x = otu.t, factors = group, sim.function = "vegdist",
              sim.method = dist, p.adjust.m = "BH",
              reduce = NULL, perm = nperm
            )
            pcutils::sanxian(dat.pairwise.adonis[, c("pairs", "R2", "p.value", "p.adjusted")],
              rows = NULL, nrow = Inf
            ) %>% print()
          }
        }
      }
    }
    if (method == "mantel") {
      # only numeric variables can do with mantel
      env %>% dplyr::select_if(\(x)is.numeric(x) & !is.factor(x)) -> env
      species.distance <- vegan::vegdist(otu.t, method = dist)
      for (i in 1:ncol(env)) {
        dd <- vegan::mantel(species.distance, vegan::vegdist(env[, i, drop = TRUE], method = dist),
          method = "pearson", permutations = nperm, na.rm = TRUE
        )
        soil <- rbind(soil, c(colnames(env)[i], dd$statistic, dd$signif))
      }
    }
    # only group variables can do with mrpp/anosim
    if (method == "anosim") {
      env %>% dplyr::select_if(\(x)class(x) == "Date" | is.factor(x) | is.character(x)) -> env
      env %>% dplyr::mutate_all(\(x)as.factor(x)) -> env
      for (i in 1:ncol(env)) {
        env[, i, drop = TRUE] -> group
        anosim_res <- vegan::anosim(otu.t, group, permutations = nperm, distance = dist)
        soil <- rbind(soil, c(colnames(env)[i], anosim_res$statistic, anosim_res$signif))
      }
    }
    if (method == "mrpp") {
      env %>% dplyr::select_if(\(x)class(x) == "Date" | is.factor(x) | is.character(x)) -> env
      env %>% dplyr::mutate_all(\(x)as.factor(x)) -> env
      for (i in 1:ncol(env)) {
        env[, i, drop = TRUE] -> group
        mrpp_res <- vegan::mrpp(otu.t, group, permutations = nperm, distance = dist)
        soil <- rbind(soil, c(colnames(env)[i], mrpp_res$A, mrpp_res$Pvalue))
      }
    }
    soil <- data.frame(group = soil[, 1], r2 = soil[, 2], p_value = soil[, 3])
  } else {
    message("Use method='adonis' to test all factors")
    message("Permanova test for all factors, notes that the order of factors would affect the result.")
    dat.div <- vegan::adonis2(otu.t ~ ., data = env, permutations = nperm, method = dist)
    dat.div <- dat.div[seq_len(nrow(dat.div) - 2), , drop = FALSE] %>% as.data.frame()
    soil <- data.frame(group = rownames(dat.div), r2 = dat.div$R2, p_value = dat.div$`Pr(>F)`)
  }
  soil$r2 <- round(as.numeric(soil$r2), 4)
  soil$p_value <- round(as.numeric(soil$p_value), 4)
  soil$sig <- soil$p_value < 0.05
  if (method %in% c("anosim", "mrpp", "mantel")) colnames(soil)[2] <- "r"
  attributes(soil)$method <- method
  class(soil) <- c("g_test", "data.frame")
  return(soil)
}

#' Multi-table test with env
#'
#' @param g_otutab multi-otutabs with first column is group
#' @param env environmental factors
#'
#' @return a mant_g object
#' @export
#'
#' @examples
#' if (requireNamespace("linkET")) {
#'   data(otutab, package = "pcutils")
#'   cbind(group = rep(c("a", "b", "c"), c(200, 100, 185)), otutab) -> g_otutab
#'   metadata[, 3:8, drop = FALSE] -> env
#'   m_group_env(g_otutab, env) -> mant_g
#'   plot(mant_g)
#' }
m_group_env <- function(g_otutab, env) {
  group <- r <- p_value <- NULL
  groups <- g_otutab$group %>% unique()
  all <- data.frame()
  for (i in groups) {
    filter(g_otutab, group == i)[, -1] -> tmp
    suppressWarnings(permanova(tmp, env, method = "mantel") -> res)
    all <- rbind(all, data.frame(spec = i, res))
  }
  all %>% dplyr::mutate(
    rd = cut(r,
      breaks = c(-Inf, 0.2, 0.4, Inf),
      labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")
    ),
    pd = cut(p_value,
      breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")
    )
  ) -> all
  all <- list(mantel_test = all, env = env)
  class(all) <- c("mant_g", class(all))
  return(all)
}

#' Plot mant_g object
#'
#' @param x mant_g object
#' @param ... add
#'
#' @return a ggplot
#' @exportS3Method
#' @method plot mant_g
#'
#' @rdname m_group_env
plot.mant_g <- function(x, ...) {
  env <- x[["env"]]
  mantel_test <- x[["mantel_test"]]
  pd <- rd <- NULL
  # if (FALSE) {
  #     lib_ps("ggcor", library = FALSE)
  #     if (isNamespaceLoaded("linkET")) lapply(c("ggcor", "linkET"), unloadNamespace)
  #     ggcor::set_scale(rev(pcutils::get_cols(pal = "bluered")),
  #         type = "gradient2n"
  #     )
  #     corp <- ggcor::quickcor(env, type = "lower") +
  #         # geom_square() +
  #         ggcor::geom_square(data = ggcor::get_data(show.diag = FALSE)) +
  #         ggcor::anno_link(aes(colour = pd, size = rd), data = mantel_test, curvature = -0.25) +
  #         # ggcor::geom_diag_label(size=3)+
  #         scale_size_manual(values = c(0.5, 1, 2)) +
  #         scale_colour_manual(values = c("#F58058", "#F7E874", "#D6D6D6")) +
  #         guides(
  #             size = guide_legend(
  #                 title = "Mantel's r",
  #                 override.aes = list(colour = "grey35"),
  #                 order = 2
  #             ),
  #             colour = guide_legend(
  #                 title = "Mantel's p",
  #                 override.aes = list(size = 3),
  #                 order = 1
  #             ),
  #             fill = guide_colorbar(title = "Pearson's r", order = 3)
  #         )
  # }
  {
    lib_ps("linkET", library = FALSE)
    corp <- linkET::qcorrplot(linkET::correlate(env), type = "lower", diag = FALSE) +
      linkET::geom_square() +
      linkET::geom_couple(aes(colour = pd, size = rd),
        data = mantel_test,
        curvature = linkET::nice_curvature()
      ) +
      scale_fill_gradientn(colours = (pcutils::get_cols(pal = "bluered"))) +
      scale_size_manual(values = c(0.5, 1, 2)) +
      scale_colour_manual(values = c("#D95F02", "#1B9E77", "#CCCCCC99")) +
      guides(
        size = guide_legend(
          title = "Mantel's r",
          override.aes = list(colour = "grey35"),
          order = 2
        ),
        colour = guide_legend(
          title = "Mantel's p",
          override.aes = list(size = 3),
          order = 1
        ),
        fill = guide_colorbar(title = "Pearson's r", order = 3)
      )
  }
  return(corp)
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
#' @return html widget
#' @examples
#' \donttest{
#' if (requireNamespace("sankeyD3") && requireNamespace("tidytree")) {
#'   data(otutab, package = "pcutils")
#'   ann_tree(taxonomy[, c(1, 5, 6, 7)], otutab) -> tree
#'   sangji_plot(tree)
#' }
#' }
sangji_plot <- function(tree, top_N = 5, notshow = c(), intermediate = FALSE, width = 3000, height = 500, ...) {
  node <- label <- parent_label <- level <- abundance <- NULL
  lib_ps("sankeyD3", library = FALSE)
  lib_ps("tidytree", library = FALSE)

  # 桑基图
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

  # le <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # le[le %in% nodes$level] -> mytax

  colnames(tree1)[(which(colnames(tree1) == "abundance") + 1):ncol(tree1)] -> mytax

  taxRank_to_depth <- setNames(seq_along(mytax) - 1, mytax)
  nodes$depth <- taxRank_to_depth[nodes$level %>% as.character()]


  links$IDsource <- match(links$parent_label, nodes$label) - 1
  links$IDtarget <- match(links$label, nodes$label) - 1
  # na.omit(links)->links

  do.call(sankeyD3::sankeyNetwork, pcutils::update_param(list(
    Links = links, Nodes = nodes,
    Source = "IDsource", Target = "IDtarget", Value = "abundance", LinkGroup = "parent_label",
    NodeID = "label", NodeGroup = "label", NodePosX = "depth", NodeValue = "abundance",
    iterations = 1000, xAxisDomain = mytax, align = "none",
    fontFamily = "arial", fontSize = 12, linkGradient = TRUE,
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
#' if (requireNamespace("plotly")) {
#'   data(otutab, package = "pcutils")
#'   ann_tree(taxonomy[, c(1, 5, 6, 7)], otutab) -> tree
#'   sunburst(tree)
#' }
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


#
#' Gene symbolid transfer to entrezIDs (human gene)
#'
#' @param genes gene symbols e.g:ASGR2
#'
#' @return gene entrezIDs dataframe
#' @export
#'
#' @examples
#' if (requireNamespace("AnnotationDbi") && requireNamespace("org.Hs.eg.db")) {
#'   genes <- c(
#'     "ASGR2", "BEST1", "SIGLEC16", "ECRP", "C1QC", "TCN2", "RNASE2",
#'     "DYSF", "C1QB", "FAM20A", "FCGR1A", "CR1", "HP", "VSIG4", "EGR1"
#'   )
#'   gene2id(genes) -> geneid
#' }
gene2id <- function(genes) {
  lib_ps("AnnotationDbi", "org.Hs.eg.db", library = FALSE)
  # entrezIDs=AnnotationDbi::mget(genes,org.Hs.egSYMBOL2EG, ifnotfound=NA)  # 找出基因对应的ID
  suppressMessages({
    entrezIDs <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID")
  })
  entrezIDs <- as.character(entrezIDs) # 获取数据
  rt <- data.frame(genes, entrezID = entrezIDs) # 添加一列entrezID
  # rt=rt[(rt[,"entrezID"])!="NA",]   # 删除没有基因的ID
  rt <- rt[!is.na(rt[, "entrezID"]), ] # 删除没有基因的ID
  return(rt)
}

#' prepare the GEO data
#'
#' @param my_id GEOid
#' @param GEO_dir GEO download dir
#' @param file the downloaded file
#'
#' @export
#' @return list(meta = meta, GSE_expr = GSE_expr)
pre_GEO <- function(my_id, GEO_dir = "GEO_data", file = NULL) {
  prode_id <- ID <- `Gene Symbol` <- GENE_SYMBOL <- Symbol <- gene_assignment <- GENE <- GeneSymbol <- n_gene <- NULL
  # 1.Importing the data
  ## change my_id to be the dataset that you want.
  lib_ps("GEOquery", "Biobase", library = FALSE)

  if (is.null(file)) {
    if (file.exists(paste0(GEO_dir, "/", my_id, "_series_matrix.txt.gz"))) {
      file <- paste0(GEO_dir, "/", my_id, "_series_matrix.txt.gz")
    }
  }

  if (is.null(file)) {
    gse <- GEOquery::getGEO(my_id[1], destdir = GEO_dir)
    gse <- gse[[1]]
  } else {
    gse <- GEOquery::getGEO(filename = file, getGPL = FALSE)
  }

  GPL_version <- gse@annotation
  if (file.exists(paste0(GEO_dir, "/", GPL_version, ".soft.gz"))) {
    GPL_data_111 <- GEOquery::getGEO(filename = paste0(GEO_dir, "/", GPL_version, ".soft.gz"))
    GPL_data_11 <- GPL_data_111@dataTable@table
  } else {
    gse <- GEOquery::getGEO(filename = file, getGPL = TRUE)
  }

  meta <- Biobase::pData(gse) ## print the sample information
  GPL_data_11 <- Biobase::fData(gse) ## print the gene annotation
  GSE_data_expr <- Biobase::exprs(gse) ## print the expression data

  # 进行注释
  GSE_data_expr <- GSE_data_expr %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::rename(prode_id = "rowname") %>%
    mutate(prode_id = as.character(prode_id))

  if (GPL_version %in% c("GPL570", "GPL571", "GPL96", "GSE22873", "GPL1261", "GPL81", "GPL8300")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, `Gene Symbol`) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "Gene Symbol")
    GPL$GeneSymbol <- gsub(" ///.*", "", GPL$GeneSymbol) # 一对多取第一个？
  } else if (GPL_version %in% c("GPL6254", "GPL10787")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, `GENE_SYMBOL`) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "GENE_SYMBOL")
    GPL$GeneSymbol <- gsub(" ///.*", "", GPL$GeneSymbol) # 一对多取第一个？
  } else if (GPL_version %in% c("GPL6947", "GPL6883", "GPL4866", "GPL4006")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, `Symbol`) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "Symbol")
    GPL$GeneSymbol <- gsub("/.*", "", GPL$GeneSymbol) # 一对多取第一个？
  } else if (GPL_version %in% c("GPL6244", "GPL17586", "GPL6246")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, gene_assignment) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "gene_assignment")
    tmp <- strsplit(GPL$GeneSymbol, split = " // ")
    lapply(tmp, \(x)x[2]) %>% do.call(c, .) -> GPL$GeneSymbol
    na.omit(GPL) -> GPL
  } else if (GPL_version %in% c("GPL1536")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, GENE) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "GENE")
    GPL$GeneSymbol <- gsub("/.*", "", GPL$GeneSymbol)
  } else if (GPL_version %in% c("GPL10739")) {
    GPL <- GPL_data_11 %>%
      dplyr::select(ID, gene_assignment) %>%
      dplyr::rename(prode_id = "ID", GeneSymbol = "gene_assignment")
    tmp <- strsplit(GPL$GeneSymbol, split = " // ")
    lapply(tmp, \(x)x[2]) %>% do.call(c, .) -> GPL$GeneSymbol
    na.omit(GPL) -> GPL
    tmp <- strsplit(GPL$GeneSymbol, split = " /// ")
    lapply(tmp, \(x)x[1]) %>% do.call(c, .) -> GPL$GeneSymbol
    na.omit(GPL) -> GPL
  } else {
    stop("unknown GPL version, please check and modify the code.")
  }

  GPL <- GPL %>%
    add_count(GeneSymbol, name = "n_gene") %>%
    add_count(prode_id, name = "n_prode_id") %>%
    dplyr::filter(n_gene < 100) %>% # 选择能有对应基因的探针,太多了说明非特异
    dplyr::select(c(-3, -4)) %>%
    mutate_all(as.character)

  GSE_data_expr[is.na(GSE_data_expr)] <- 0
  GSE_expr <- GSE_data_expr %>%
    inner_join(GPL, ., by = "prode_id") %>% # p2s和上一步得到的结果再取交集，p2s放在右边是以它为准
    dplyr::select(-1) %>% # 去除第一列probe_id
    data.frame() %>% # 因为aggregate需要数据框格式
    aggregate(. ~ GeneSymbol, data = ., mean) %>% # 以symbol作为因子水平，把相似的数据放在一起取均值，最大值max，中位值median
    tibble::column_to_rownames(var = "GeneSymbol")

  rownames(GSE_expr) %>% gsub(" ", "", .) -> rownames(GSE_expr)

  saveRDS(list(meta = meta, GSE_expr = GSE_expr), file = paste0(GEO_dir, "/", my_id, ".RDS"))
  return(list(meta = meta, GSE_expr = GSE_expr))
}
