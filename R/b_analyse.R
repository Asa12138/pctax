# b_diversity============

dist_3col <- function(dist) {
  dist <- as.matrix(dist)
  rowname <- rownames(dist)
  colname <- colnames(dist)
  rown <- row(dist)
  coln <- col(dist)
  dist.v <- as.vector(stats::as.dist(dist))
  rown.v <- as.vector(stats::as.dist(rown))
  coln.v <- as.vector(stats::as.dist(coln))
  res <- data.frame(
    name1 = rowname[rown.v], name2 = colname[coln.v],
    dis = dist.v
  )
  res
}

#' Calculate distance for otutab
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...see \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @param spe_nwk a phylo tree if use unifrac...
#'
#' @return dist
#' @export
#' @examples
#' data(otutab, package = "pcutils")
#' mat_dist(otutab)
mat_dist <- function(otutab, method = "bray", spe_nwk = NULL) {
  totu <- t(otutab)
  if (is.null(spe_nwk)) {
    vegan::vegdist(totu, method = method) -> a
  } else {
    lib_ps("picante", library = FALSE)
    picante::match.phylo.comm(spe_nwk, totu) -> match_p
    match_p$phy -> spe_nwk
    match_p$comm -> totu
    a <- switch(method,
      "unifrac" = picante::unifrac(totu, spe_nwk),
      "b_mpd" = picante::comdist(totu, stats::cophenetic(spe_nwk)),
      "b_mntd" = picante::comdistnt(totu, stats::cophenetic(spe_nwk)),
      "phylosor" = picante::phylosor(totu, spe_nwk)
      # "pcd"= picante::pcd(totu, spe_nwk),
      # picante::psd(samp,test)
      # picante::raoD(samp,test)
    )
  }
  return(a)
}

#' Transfer dist to b_dist
#'
#' @param dist a dist object
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return a b_dist with annotation by group
#' @export
#' @examples
#' data(otutab, package = "pcutils")
#' mat_dist(otutab) %>% as.b_dist(., group_df = metadata["Group"]) -> aa
#' plot(aa)
#' plot(aa, mode = 2)
as.b_dist <- function(dist, group_df = NULL) {
  group1 <- group2 <- variable <- NULL
  # 将dist矩阵转化为group注释的b_dist对象
  stopifnot(inherits(dist, "dist"))
  dist_3col(dist) -> aa
  if (!is.null(group_df)) {
    stopifnot(is.data.frame(group_df))
    group <- group_df %>% unlist()
    if (FALSE) {
      # 另一种方法
      as.matrix(a) -> b
      group <- droplevels.factor(group)
      utils::combn(levels(factor(group)), 2) %>% apply(., 2, FUN = function(x) paste0(x[1], "_", x[2])) -> pair
      list() -> b_ls
      for (i in 1:nrow(b)) {
        for (j in 1:nrow(b)) {
          b_ls[[paste0(group[i], "_", group[j])]] <- c(b_ls[[paste0(group[i], "_", group[j])]], b[i, j])
        }
      }
      b_ls <- b_ls[pair]
      suppressMessages(as.data.frame(b_ls) %>% melt() -> aa)
      colnames(aa)[2] <- "dis"
    }
    com <- \(group1, group2, levels){
      factor(c(group1, group2), levels = levels) %>%
        sort() %>%
        paste(., collapse = "_")
    }
    if (nrow(as.matrix(dist)) != length(group)) stop("group length wrong")
    aa <- dplyr::left_join(aa, data.frame(name1 = rownames(group_df), group1 = group), by = "name1")
    aa <- dplyr::left_join(aa, data.frame(name2 = rownames(group_df), group2 = group), by = "name2")
    aa %>% dplyr::mutate(a = com(group1, group2, levels(group)))

    aa$variable <- apply(aa, 1, function(x) com(x[4], x[5], levels(group)))
    aa %>% dplyr::mutate(group = ifelse(group1 == group2, "intra", "inter")) -> aa
    aa %>% dplyr::mutate(variable = ifelse(group == "inter", variable, as.character(group1))) -> aa
  }
  class(aa) <- c("b_dist", class(aa))
  return(aa)
}

#' Transfer b_dist to dist
#'
#' @param m a b_dist object
#' @param diag 	logical value indicating whether the diagonal of the distance matrix should be printed by \code{print.dist}.
#' @param upper logical value indicating whether the upper triangle of the distance matrix should be printed by \code{print.dist}.
#'
#' @return dist
#' @importFrom stats as.dist
#' @exportS3Method
#' @method as.dist b_dist
as.dist.b_dist <- function(m, diag = FALSE, upper = FALSE) {
  reshape2::dcast(m, name1 ~ name2, value.var = "dis") %>% column_to_rownames("name1") -> tmp
  rbind(rep(NA, ncol(tmp)), tmp) -> tmp
  rownames(tmp)[1] <- setdiff(colnames(tmp), rownames(tmp))
  cbind(tmp, rep(NA, nrow(tmp))) -> tmp
  colnames(tmp)[ncol(tmp)] <- setdiff(rownames(tmp), colnames(tmp))
  tmp[rownames(tmp), rownames(tmp)] -> tmp
  tmp[is.na(tmp)] <- 0
  tmp + t(tmp) -> tmp
  return(stats::as.dist(tmp, diag = diag, upper = upper))
}

#' Plot dist
#'
#' @param x a dist
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param ... additional
#'
#' @return a pheatmap
#' @exportS3Method
#' @method plot dist
#' @rdname as.b_dist
plot.dist <- function(x, group_df = NULL, ...) {
  lib_ps("pheatmap", library = FALSE)
  ab <- as.matrix(x) %>% as.data.frame()
  if (!is.null(group_df)) {
    group <- group_df %>% unlist()
    ann <- group_df %>% dplyr::arrange(factor(group))
    ab <- ab[rownames(ann), rownames(ann)]
    pheatmap::pheatmap(ab,
      annotation_row = ann, annotation_col = ann, cluster_rows = FALSE,
      cluster_cols = FALSE, border_color = FALSE, ...
    )
  } else {
    pheatmap::pheatmap(ab)
  }
}

#' Plot b_dist
#'
#' @param x a b_dist
#' @param mode 1~3
#' @param c_group "inter" or "intra" or both to plot
#' @param ... additional
#'
#' @import ggplot2
#' @return a ggplot or pheatmap
#' @exportS3Method
#' @rdname as.b_dist
#'
plot.b_dist <- function(x, mode = 1, c_group = "inter", ...) {
  aa <- x
  group <- dis <- NULL
  if (mode == 1) {
    aa$variable <- factor(aa$variable, levels = unique(c(levels(aa$group1), unique(aa$variable))))
    dplyr::filter(aa, group %in% c_group) -> aa1
    p <- pcutils::group_box(aa1$dis, group = aa1$variable, ...)
    return(p)
  }
  if (mode == 2) {
    pl <- ggplot(aa, aes(x = dis, col = group, fill = group)) +
      geom_density(alpha = 0.5) +
      geom_rug() +
      geom_vline(aes(xintercept = median, col = group),
        data = aa %>% dplyr::group_by(group) %>%
          dplyr::summarise(median = median(dis)), linetype = 2
      ) +
      scale_color_manual(values = c("#E88493", "#51A4E0")) +
      scale_fill_manual(values = c("#E88493", "#51A4E0")) +
      labs(x = "Bray-Curtis Distance")
    # 获得ggplot绝对画板大小
    (stats::wilcox.test(aa$dis ~ aa$group) -> wt)
    ggplot_build(pl) -> plot
    p <- pl + annotate("text",
      x = 0.9 * max(plot$data[[1]]$x), y = max(plot$data[[1]]$y),
      label = paste0("p=", signif(wt$p.value, 4))
    )
    return(p)
  }
  if (mode == 3) {
    ann <- data.frame(name = c(aa$name1, aa$name2), group = c(aa$group1, aa$group2)) %>%
      distinct() %>%
      column_to_rownames("name")
    plot.dist(as.dist.b_dist(aa), ann)
  }
}


#' Lm for sample similarity and geographical distance
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param geo a two-columns dataframe, first is latitude, second is longitude
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...see \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @param spe_nwk a phylo tree if use unifrac...
#' @param ... additional
#'
#' @import ggplot2
#' @return a ggplot
#' @export
#' @references
#' Graco-Roza, C. et al. (2022) Distance decay 2.0 - A global synthesis of taxonomic and functional turnover in ecological communities. Glob Ecol Biogeogr 31, 1399–1421.
#' @examples
#' if (requireNamespace("NST") && requireNamespace("geosphere")) {
#'   library(ggplot2)
#'   data(otutab, package = "pcutils")
#'   metadata[, c("lat", "long")] -> geo
#'   geo_sim(otutab, geo) -> geo_res
#'   pcutils::my_lm(geo_res[4], "dis.geo", geo_res)
#' }
geo_sim <- function(otutab, geo, method = "bray", spe_nwk = NULL, ...) {
  lib_ps("NST", "geosphere", library = FALSE)
  dis <- NULL
  # 经纬度数据转换
  # 直接欧式距离算不太准
  # pcutils::toXY(geo)%>%ggscatter(.,"X","Y",label = rownames(geo))
  geosphere::distm(geo[, 2:1]) -> geo_dist
  rownames(geo_dist) <- colnames(geo_dist) <- rownames(geo)

  geo_dist %>%
    as.dist() %>%
    dist_3col() %>%
    dplyr::mutate(dis = dis / 1000) -> geo_dist

  # 这一步可以选择各种b_diversity指数，甚至是谱系与功能多样性指数
  mat_dist(otutab, method, spe_nwk) %>%
    dist_3col() %>%
    dplyr::mutate(dis = 1 - dis) -> similarity

  merge(geo_dist, similarity, by = c("name1", "name2"), suffixes = c(".geo", ".b")) -> a
  return(a)
}

#' Beta_diversity Ordination: dimensionality reduction
#'
#' @description
#' Species abundance data can be preprocessed with Hellinger transformation or chord transformation data before PCA analysis. Because the Hellinger distance or chord distance with-without data is equal to \eqn{\sqrt2\sqrt{1-Ochiai\ similarity}}, therefore, the sorting diagram (type 1 scale) of PCA analysis after Hellinger transformation or chord transformation with-without data is internal sample The distance between the squares is the Ochiai distance. \eqn{\sqrt2\sqrt{1-Ochiai\ similarity}} is a distance measure, which is also suitable for the analysis of species data. The processed data is then used for pca without norm.
#'
#' @param otutab object
#' @param ... additional
#'
#' @return b_res object
#' @export
#' @references
#' <https://www.jianshu.com/p/9694c0b6302d>
#' <https://zhuanlan.zhihu.com/p/25501130>
#' @examples
#' data(otutab, package = "pcutils")
#' b_analyse(otutab, method = "pca") -> b_res
#' plot(b_res, "Group", metadata)
b_analyse <- function(otutab, ...) {
  UseMethod("b_analyse")
}

#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param norm should normalized or not? (hellinger)
#' @param method one of "pca","pcoa","ca","dca","nmds","plsda","tsne","umap","lda","all"
#' @param group if needed, give a group vector
#' @param dist if use pcoa or nmds, your can choose a dist method (default: bray) or input a distance matrix.
#' @param ndim how many dimension be kept? (default:2). 3 for b_res_3d()
#' @param scale scale, default: FALSE
#' @param ... add
#'
#' @export
#' @method b_analyse data.frame
#' @rdname b_analyse
b_analyse.data.frame <- function(otutab, norm = TRUE, method = c("pca", "nmds"), group = NULL,
                                 dist = "bray", ndim = 2, scale = FALSE, ...) {
  Comp1 <- Comp2 <- CS1 <- CS2 <- NMDS1 <- NMDS2 <- PLS_DA1 <- PLS_DA2 <- NULL
  lib_ps("ade4", library = FALSE)
  all <- c("pca", "pcoa", "ca", "dca", "nmds", "plsda", "tsne", "umap", "lda", "all")
  if (!all(method %in% all)) stop("method is one of ", paste(all, collapse = ","))
  if ("all" %in% method) method <- all[-length(all)]

  data.frame(t(otutab), check.names = FALSE) -> dat
  if (norm) {
    dat.h <- vegan::decostand(dat, "hellinger", MARGIN = 1)
  } else {
    dat.h <- dat
  }

  # storage dataframes
  data.frame(name = colnames(otutab)) -> sample_site
  data.frame(eigs = paste("eig", 1:ndim)) -> sample_eig
  data.frame(var = rownames(otutab)) -> var_site
  data.frame(var = rownames(otutab)) -> var_contri
  # PCA:principal components analysis 主成分分析
  if ("pca" %in% method) {
    # dat.pca <- rda(dat.h, scale = FALSE)#scale=TRUE,对列变量进行标准化
    # summary(dat.pca, scaling = 1)->a#scaling=1展示样本距离，=2展示变量夹角
    # cbind(sample_eig,a[["cont"]][["importance"]][2,1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'PCA'
    # cbind(sample_site,a$sites[,1:2])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('PCA1', 'PCA2')

    dat.pca <- ade4::dudi.pca(dat.h, scale = scale, scannf = FALSE, nf = 3)
    cbind(sample_eig, (dat.pca$eig / sum(dat.pca$eig))[1:ndim]) -> sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- "PCA"
    cbind(sample_site, dat.pca$li[, 1:ndim]) -> sample_site # 绘图坐标
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("PCA", 1:ndim)
    cbind(var_site, dat.pca$c1[, 1:2]) -> var_site # 变量坐标
    colnames(var_site)[(ncol(var_site) - 1):ncol(var_site)] <- paste0("PCA", 1:2)
    dat.pca$co[, 1:2] %>%
      dplyr::transmute(contri = ((Comp1^2 / (dat.pca$eig)[1]) + (Comp2^2 / (dat.pca$eig)[2])) * 100 / 2) %>%
      cbind(var_contri, .) -> var_contri
    colnames(var_contri)[ncol(var_contri)] <- "PCA"

    # 这个结果可以用fviz_pca直接可视化
    # if(FALSE){
    #   fviz_screeplot(dat.pca)
    #   fviz_pca_ind(dat.pca, col.ind = "cos2")
    #   fviz_pca_var(dat.pca)
    # }
  }
  # PcoA:principal coordinate analysis 主坐标分析
  if ("pcoa" %in% method) {
    # dat.pcoa <- ape::pcoa(vegdist(dat.h,method = dist))
    # cbind(sample_eig,dat.pcoa$values$Relative_eig[1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'PCoA'
    # cbind(sample_site,dat.pcoa$vectors[,c(1,2)])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('PCoA1', 'PCoA2')

    if (is.matrix(dist) | is.data.frame(dist)) {
      message("use the `dist` as a distance matrix.")
      dist_df <- as.dist(dist)
    } else if (inherits(dist, "dist")) {
      message("use the `dist` as a distance.")
      dist_df <- dist
    } else {
      dist_df <- vegan::vegdist(dat.h, method = dist, na.rm = TRUE)
    }

    dat.pco <- ade4::dudi.pco(dist_df, scannf = FALSE, nf = 3)
    cbind(sample_eig, (dat.pco$eig / sum(dat.pco$eig))[1:ndim]) -> sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- "PCoA"
    cbind(sample_site, dat.pco$li[, 1:ndim]) -> sample_site # 绘图坐标
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("PCoA", 1:ndim)
    # cbind(var_site,dat.pco$c1[,1:2])->var_site#变量坐标
    # colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('PCoA1', 'PCoA2')
    # dat.pco$co[,1:2]%>%transmute(contri=((Comp1^2/(dat.pco$eig)[1])+(Comp2^2/(dat.pco$eig)[2]))*100/2)%>%
    #   cbind(var_contri,.)->var_contri
    # colnames(var_contri)[ncol(var_contri)] <- 'PCoA'
  }
  # CA:correspondence analysis 对应分析
  if ("ca" %in% method) {
    # dat.ca <- cca(dat.h,scale = FALSE)
    # summary(dat.ca, scaling = 1)->b
    # cbind(sample_eig,(dat.ca$CA$eig/sum(dat.ca$CA$eig))[1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'CA'
    # cbind(sample_site,b$sites[,1:2])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('CA1', 'CA2')

    dat.coa <- ade4::dudi.coa(dat.h, scannf = FALSE, nf = 3)
    cbind(sample_eig, (dat.coa$eig / sum(dat.coa$eig))[1:ndim]) -> sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- "CA"
    cbind(sample_site, dat.coa$li[, 1:ndim]) -> sample_site # 绘图坐标
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("CA", 1:ndim)
    cbind(var_site, dat.coa$co[, 1:2]) -> var_site # 变量坐标
    colnames(var_site)[(ncol(var_site) - 1):ncol(var_site)] <- c("CA1", "CA2")
    dat.coa$c1[, 1:2] %>%
      dplyr::transmute(contri = ((CS1^2 / (dat.coa$eig)[1]) + (CS2^2 / (dat.coa$eig)[2])) * 100 / 2) %>%
      cbind(var_contri, .) -> var_contri
    colnames(var_contri)[ncol(var_contri)] <- "CA"
  }
  # DCA:detrended correspondence analysis 去趋势分析
  if ("dca" %in% method) {
    dat.dca <- vegan::decorana(dat.h)
    # all.dca=summary(dat.dca)
    cbind(sample_eig, (dat.dca$evals / sum(dat.dca$evals))[1:ndim]) -> sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- "DCA"
    cbind(sample_site, dat.dca$rproj[, 1:ndim]) -> sample_site # 绘图坐标
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("DCA", 1:ndim)
    # cbind(var_site,dat.pco$c1[,1:2])->var_site#变量坐标
    # colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('DCA1', 'DCA2')
    # dat.pco$co[,1:2]%>%transmute(contri=((Comp1^2/(dat.pco$eig)[1])+(Comp2^2/(dat.pco$eig)[2]))*100/2)%>%
    #   cbind(var_contri,.)->var_contri
    # colnames(var_contri)[ncol(var_contri)] <- 'DCA'
  }
  # NMDS:non-metric multi-dimensinal scaling 非度量多维尺度分析
  if ("nmds" %in% method) {
    suppressMessages(dat.nmds <- vegan::metaMDS(dat.h, k = ndim, distance = dist))
    cbind(sample_eig, data.frame(a = rep("", ndim))) -> sample_eig # 两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- "NMDS"
    # cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site, dat.nmds$points) -> sample_site
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("NMDS", 1:ndim)
    cbind(var_site, dat.nmds$species[, 1:2]) -> var_site # 变量坐标
    colnames(var_site)[(ncol(var_site) - 1):ncol(var_site)] <- c("NMDS1", "NMDS2")
    cbind(var_contri, var_site %>% transmute(((NMDS1^2 + NMDS2^2) / 2))) -> var_contri
    colnames(var_contri)[ncol(var_contri)] <- "NMDS"
  }
  # PLS-DA:partial least squares discriminant analysis 偏最小二乘法判别分析
  if ("plsda" %in% method) {
    if (is.null(group)) stop("plsda need group!")
    lib_ps("mixOmics", library = FALSE)
    mixOmics::plsda(dat.h, group, ncomp = ndim) -> plsdat
    cbind(sample_eig, data.frame(plsdat = plsdat$prop_expl_var$X)) -> sample_eig # 两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- "PLS_DA"
    # cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site, plsdat$variates$X) -> sample_site
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("PLS_DA", 1:ndim)
    cbind(var_site, plsdat$mat.c[, 1:2] * 100) -> var_site # 变量坐标
    colnames(var_site)[(ncol(var_site) - 1):ncol(var_site)] <- c("PLS_DA1", "PLS_DA2")
    var_site[, (ncol(var_site) - 1):ncol(var_site)] %>%
      dplyr::transmute(contri = ((PLS_DA1^2 / (plsdat$prop_expl_var$X)[1]) + (PLS_DA2^2 / (plsdat$prop_expl_var$X)[2])) * 100 / 2) %>%
      cbind(var_contri, .) -> var_contri
    colnames(var_contri)[ncol(var_contri)] <- "PLS_DA"
  }
  # LDA:linear discriminant analysis 线性判别分析
  if ("lda" %in% method) {
    if (is.null(group)) stop("lda need group!")
    lib_ps("MASS", library = FALSE)
    lda.sol <- MASS::lda(dat.h, group)
    P <- lda.sol$scaling
    # 将均值向量降维
    means <- lda.sol$means %*% P
    # 加权平均的出总的降维均值向量，权重就是lda.sol$prior
    total_means <- as.vector(lda.sol$prior %*% means)
    # 把样本降维并平移
    x <- as.matrix(dat.h) %*% P - (rep(1, nrow(dat.h)) %o% total_means)
    # ldat=(lda.sol$svd)**2/sum((lda.sol$svd)**2)#两轴的解释度
    cbind(sample_eig, a = "") -> sample_eig # 两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- "LDA"
    # cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site, x) -> sample_site
    colnames(sample_site)[(ncol(sample_site) - 1):ncol(sample_site)] <- c("LDA1", "LDA2")
  }
  # t-SNE:t-distributed stochastic neighbor embedding t-分布式随机邻域嵌入
  if ("tsne" %in% method) {
    lib_ps("Rtsne", library = FALSE)
    Rtsne::Rtsne(dat.h,
      dims = ndim, pca = TRUE,
      max_iter = 1000,
      theta = 0.4,
      perplexity = 5, # 默认20，会报错
      verbose = FALSE
    ) -> tsne
    cbind(sample_eig, a = "") -> sample_eig # 两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- "t_SNE"
    # cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site, tsne$Y) -> sample_site
    colnames(sample_site)[(ncol(sample_site) - ndim + 1):ncol(sample_site)] <- paste0("t_SNE", 1:ndim)
  }
  # UMAP:uniform manifold approximation and projection 均一流形近似投影
  if ("umap" %in% method) {
    lib_ps("umap", library = FALSE)
    umap::umap(dat.h) -> umapr
    cbind(sample_eig, a = "") -> sample_eig # 两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- "UMAP"
    # cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site, umapr$layout) -> sample_site
    colnames(sample_site)[(ncol(sample_site) - 1):ncol(sample_site)] <- c("UMAP1", "UMAP2")
  }
  # # fso: fuzzy set ordination 模糊集排序
  # if (FALSE) {
  #     lib_ps("fso", library = FALSE)
  #     env.fso <- fso::fso(metadata[, 12:13], vegdist(dat.h))
  #     data.frame(x = env.fso$data, y = env.fso$mu) %>% ggplot(., aes(x = y.1, y = y.2)) +
  #         geom_point(aes(col = metadata$Group))
  # }
  # # sofm:self-organizing feature map 自组织特质
  # if (FALSE) {
  #     lib_ps("kohonen", library = FALSE)
  #     # kohonen::xyf()
  # }

  message("four dataframes in a list, 1 is eig, 2 is sample_site, 3 is var, 4 is var contribution")
  b_res <- list(sample_eig = sample_eig, sample_site = sample_site, var_site = var_site, var_contri = var_contri)
  class(b_res) <- "b_res"
  return(b_res)
}

is.continuous <- function(x) {
  is.numeric(x) | inherits(x, "Date")
}


# base plot for pca/rda
plot_b_like <- function(plotdat, mode = 1, pal = NULL, sample_label = TRUE, stat_ellipse = TRUE, groupname = "level", groupname2 = "level2", ...) {
  x1 <- x2 <- level <- level2 <- x1.cen <- x2.cen <- NULL
  if (mode == 1) {
    plist <- {
      ggplot(plotdat, aes(x = x1, y = x2)) +
        do.call(
          geom_point,
          update_param(list(mapping = aes(fill = level, shape = level2), color = "black", size = 2), list(...))
        ) +
        # 可在这里修改点的透明度、大小
        scale_shape_manual(values = 21:25) +
        guides(fill = guide_legend(override.aes = list(shape = 21))) +
        geom_vline(xintercept = 0, color = "gray", linewidth = 0.4) +
        geom_hline(yintercept = 0, color = "gray", linewidth = 0.4)
    }

    if (!is.continuous(plotdat$level) & (stat_ellipse == 1)) {
      plist <- plist + stat_ellipse(aes(fill = level), type = "norm", geom = "polygon", alpha = 0.1, color = NA, level = 0.68)
    }
  } else if (mode == 2) {
    plist <- {
      ggplot(plotdat, aes(x = x1, y = x2)) +
        do.call(
          geom_point,
          update_param(list(mapping = aes(color = level, shape = level2), size = 2), list(...))
        ) +
        # 可在这里修改点的透明度、大小
        geom_vline(xintercept = 0, color = "gray", linewidth = 0.4) +
        geom_hline(yintercept = 0, color = "gray", linewidth = 0.4)
    }
    if (!is.continuous(plotdat$level) & (stat_ellipse == 1)) {
      plist <- plist + stat_ellipse(aes(color = level), level = 0.68)
    }
  } else if (mode == 3) {
    if (is.continuous(plotdat$level)) stop("Group is continous! try mode 1 or mode 2")
    centroid <- aggregate(cbind(x1, x2) ~ level,
      data = plotdat, FUN = mean
    )
    plotdat1 <- dplyr::left_join(plotdat, centroid, by = "level", suffix = c("", ".cen"))
    plotdat <- data.frame(plotdat1, row.names = rownames(plotdat))
    plist <- {
      ggplot(plotdat, aes(x = x1, y = x2)) +
        do.call(
          geom_point,
          update_param(list(mapping = aes(color = level, shape = level2), size = 2), list(...))
        ) +
        # 可在这里修改点的透明度、大小
        geom_vline(xintercept = 0, color = "gray", linewidth = 0.4) +
        geom_hline(yintercept = 0, color = "gray", linewidth = 0.4) +
        geom_segment(aes(xend = x1.cen, yend = x2.cen, color = level), show.legend = FALSE) +
        geom_label(
          data = centroid, aes(label = level, fill = level), size = 5,
          show.legend = FALSE, color = "white"
        )
    }
  }
  # sample_label
  if (sample_label) {
    lib_ps("ggrepel", library = FALSE)
    plist <- plist + ggrepel::geom_text_repel(aes(x = x1, y = x2, label = rownames(plotdat)),
      col = "black", size = 2.5
    )
  }

  if (!is.continuous(plotdat$level)) {
    plist <- plist +
      scale_color_manual(values = pal, name = groupname) + # 可在这里修改点的颜色
      scale_fill_manual(values = pal, name = groupname)
  } else if (inherits(plotdat$level, "Date")) {
    if (mode == 1) {
      plist <- plist +
        scale_fill_gradientn(colours = pal, trans = scales::date_trans(), name = groupname)
    } else {
      plist <- plist +
        scale_color_gradientn(colours = pal, trans = scales::date_trans(), name = groupname)
    }
  } else {
    plist <- plist +
      scale_color_gradientn(colours = pal, name = groupname) +
      scale_fill_gradientn(colours = pal, name = groupname)
  }

  if (groupname2 == "_group2") {
    plist <- plist + guides(shape = guide_none())
  } else {
    plist <- plist + guides(shape = guide_legend(title = groupname2))
  }

  if (!is.continuous(plotdat$level2) & (stat_ellipse == 2)) {
    if (mode == 1) {
      plist <- plist +
        ggnewscale::new_scale_fill() +
        stat_ellipse(aes(fill = level2), type = "norm", geom = "polygon", alpha = 0.1, color = NA, level = 0.68) +
        scale_fill_manual(name = groupname2, values = pcutils::get_cols(nlevels(factor(plotdat$level2))))
    }
    if (mode == 2) {
      plist <- plist +
        ggnewscale::new_scale_color() +
        stat_ellipse(aes(color = level2), level = 0.68) +
        scale_color_manual(name = groupname2, values = pcutils::get_cols(nlevels(factor(plotdat$level2))))
    }
  }
  plist <- plist + pctax_theme
  return(plist)
}

#' Plot a b_res
#'
#' @param x a b_res object
#' @param Group group vector for color
#' @param metadata metadata contain Group
#' @param Group2 mapping point shape
#' @param mode plot mode:1~3
#' @param bi plot variables segments?
#' @param Topn how many variables to show?
#' @param rate segments length rate
#' @param margin plot the margin boxplot?
#' @param box_margin margin plot box or density?
#' @param pal colors for group
#' @param sample_label plot the labels of samples?
#' @param coord_fix fix the coordinates y/x ratio
#' @param bi_text_size biplot text size
#' @param ... add
#' @param stat_ellipse plot the stat_ellipse?
#' @param box_param box_param for \code{\link[pcutils]{group_box}}
#' @param margin_label margin_label, TRUE
#' @param permanova_res permanova result
#' @param text_param text_param for \code{\link[ggplot2]{annotate}}
#'
#' @return a ggplot
#' @exportS3Method
#' @method plot b_res
#'
#' @seealso \code{\link{b_analyse}}
#'
plot.b_res <- function(x, Group, metadata = NULL, Group2 = NULL,
                       mode = 1, bi = FALSE, Topn = 10, rate = 1, margin = FALSE,
                       margin_label = TRUE, permanova_res = NULL, text_param = list(),
                       box_margin = TRUE, box_param = list(), pal = NULL, sample_label = TRUE,
                       stat_ellipse = TRUE, coord_fix = FALSE,
                       bi_text_size = 3, ...) {
  b_res <- x
  contri <- x1 <- x2 <- level <- NULL
  # prepare metadata
  if (!is.null(metadata)) {
    if (!Group %in% colnames(metadata)) stop("Group should be one of colnames(metadata)")
    if (!is.null(Group2)) {
      if (!(Group2 %in% colnames(metadata))) stop("Group2 should be one of colnames(metadata)")
    }
    idx <- rownames(metadata) %in% rownames(b_res$sample_site)
    if (!all(rownames(metadata) %in% rownames(b_res$sample_site))) message("rownames don't match in b_res and metadata")
    metadata <- metadata[idx, , drop = FALSE]
    b_res$sample_site <- b_res$sample_site[rownames(metadata), , drop = FALSE]
  } else {
    metadata <- data.frame(row.names = rownames(b_res$sample_site), group = Group)
    Group <- "group"
    # input group2?
    if (is.null(Group2)) {
      Group2 <- "_shape"
      metadata <- data.frame(metadata, `_group2` = Group2, check.names = FALSE)
      Group2 <- "_group2"
    } else {
      metadata <- data.frame(metadata, group2 = Group2)
      Group2 <- "group2"
    }
  }
  # input group2?
  if (is.null(Group2)) {
    Group2 <- "_shape"
    metadata <- data.frame(metadata, `_group2` = Group2, check.names = FALSE)
    Group2 <- "_group2"
  }

  if (is.null(pal)) {
    if (is.continuous(metadata[, Group])) {
      pal <- RColorBrewer::brewer.pal(8, "Reds")
    } else {
      pal <- pcutils::get_cols(n = length(unique(metadata[, Group])), pal = RColorBrewer::brewer.pal(5, "Set2"))
    }
  }

  # mode 代表用哪种风格画图，常用的1-3已经准备好了，临时改的话添加4就行。
  plist <- list()
  for (i in colnames(b_res$sample_eig)[-1]) {
    b_res$sample_site %>% dplyr::select(dplyr::starts_with(i)) -> tmp
    plotdat <- data.frame(x1 = tmp[, 1], x2 = tmp[, 2], level = metadata[, Group], level2 = metadata[, Group2], row.names = b_res$sample_site$name)
    b_res$sample_eig %>%
      dplyr::select(dplyr::starts_with(i)) %>%
      unlist() -> eig

    plist[[i]] <- plot_b_like(plotdat, mode = mode, pal = pal, sample_label = sample_label, stat_ellipse = stat_ellipse, groupname = Group, groupname2 = Group2, ...)

    # add permanova result
    if (!is.null(permanova_res)) {
      if (inherits(permanova_res, "g_test")) {
        p_res <- permanova_res[which(permanova_res$group == Group), ]
        if (attributes(p_res)$method == "adonis") {
          p_label <- paste0("Adonis\nR2=", round(p_res$r2, 4), "\nP=", round(p_res$p_value, 4))
        } else if (attributes(p_res)$method == "anosim") {
          p_label <- paste0("Anosim\nR=", round(p_res$r, 4), "\nP=", round(p_res$p_value, 4))
        } else if (attributes(p_res)$method == "mrpp") {
          p_label <- paste0("MRPP\nR=", round(p_res$r, 4), "\nP=", round(p_res$p_value, 4))
        } else if (attributes(p_res)$method == "mantel") {
          p_label <- paste0("Mantel\nR=", round(p_res$r, 4), "\nP=", round(p_res$p_value, 4))
        } else if (attributes(p_res)$method == "envfit") {
          p_label <- paste0("Envfit\nR2=", round(p_res$r2, 4), "\nP=", round(p_res$p_value, 4))
        } else {
          p_label <- paste0(attributes(p_res)$method, "\nP=", round(p_res$p_value, 4))
        }

        # q: 判断ggplot版本大于3.5.0
        if (utils::packageVersion("ggplot2") >= "3.5.0") {
          plist[[i]] <- plist[[i]] + do.call(annotate, pcutils::update_param(
            list(
              geom = "text", x = I(0.88), y = I(0.88),
              label = p_label, size = 3, color = "black", fontface = "bold"
            ),
            text_param
          ))
        } else {
          lims <- pcutils::ggplot_lim(plist[[i]])
          lapply(lims, \(i)i[1] + 0.93 * (i[2] - i[1])) -> lims
          plist[[i]] <- plist[[i]] + do.call(annotate, pcutils::update_param(
            list(
              geom = "text", x = lims$x, y = lims$y,
              label = p_label, size = 3, color = "black", fontface = "bold"
            ),
            text_param
          ))
        }
      } else {
        warning("permanova_res should be a g_test object")
      }
    }

    if (bi == TRUE) {
      b_res$var_site %>% dplyr::select(dplyr::starts_with(i)) -> tmp
      tmp * rate -> tmp
      if ((b_res$var_site %>% dplyr::select(dplyr::starts_with(i)) %>% ncol()) > 0) {
        cbind(var = b_res$var_site[, 1], tmp, contri = b_res$var_contri[, i]) -> bi_tmp
        colnames(bi_tmp) <- c("var", "x1", "x2", "contri")
        bi_tmp %>% dplyr::top_n(Topn, contri) -> bi_tmp1
        plist[[i]] <- plist[[i]] +
          ggnewscale::new_scale_color() +
          geom_segment(
            data = bi_tmp1, aes(x = 0, y = 0, xend = x1, yend = x2, color = contri),
            size = 0.3, arrow = arrow(length = unit(0.1, "inches"))
          ) +
          ggrepel::geom_text_repel(data = bi_tmp1, aes(
            x = x1 * 1.05, y = x2 * 1.05,
            color = contri, label = var
          ), size = bi_text_size) +
          scale_color_gradientn(name = "Contribution", colours = c("#00AFBB", "#E7B800", "#FC4E07"))
      }
    }

    # labs on axis
    if (i %in% c("PCA", "PCoA", "CA", "PLS_DA")) {
      plist[[i]] <- plist[[i]] +
        labs(
          x = paste(paste(i, "1: ", sep = ""), round(100 * eig[1], 2), "%"),
          y = paste(paste(i, "2: ", sep = ""), round(100 * eig[2], 2), "%")
        )
      if (coord_fix) plist[[i]] <- plist[[i]] + coord_fixed(sqrt(eig[2] / eig[1]))
    } else {
      plist[[i]] <- plist[[i]] + labs(x = paste0(i, "1"), y = paste0(i, "2"))
      if (coord_fix) plist[[i]] <- plist[[i]] + coord_fixed(1.2)
    }

    if (margin) {
      plist[[i]] <- b_plot_margin(plist[[i]], plotdat, box_margin = box_margin, pal = pal, box_param = box_param, margin_label = margin_label)
    }
  }
  if (length(plist) == 1) plist <- plist[[1]]
  return(plist)
}


b_plot_margin <- function(p, plotdat, pal, box_margin, box_param, margin_label, ...) {
  lib_ps("aplot", library = FALSE)
  x2 <- level <- x1 <- NULL
  common_theme <- theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.line = element_blank(), axis.title = element_blank(),
      plot.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  if (box_margin) {
    p1 <- do.call(pcutils::group_box, pcutils::update_param(list(tab = plotdat["x2"], group = "level", metadata = plotdat), box_param))
    p1 <- p1 +
      common_theme +
      theme(
        axis.line.x = element_line(),
        axis.ticks.x = element_line()
      )
    p2 <- do.call(pcutils::group_box, pcutils::update_param(list(tab = plotdat["x1"], group = "level", metadata = plotdat), box_param))
    p2 <- p2 +
      common_theme +
      theme(
        axis.line.y = element_line(),
        axis.ticks.y = element_line()
      ) + coord_flip()

    if (margin_label) {
      p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = rel(0.8)))
      p2 <- p2 + theme(axis.text.y = element_text(size = rel(0.8)))
    }

    if (any(c("p_value1", "p_value2", "alpha") %in% names(box_param))) {
      if (any(unlist(box_param[intersect(c("p_value1", "p_value2", "alpha"), names(box_param))]))) {
        pcutils::ggplot_lim(p) -> main_p_lim
        pcutils::ggplot_lim(p1) -> p1_lim
        pcutils::ggplot_lim(p2) -> p2_lim
        if (p1_lim$y[2] > main_p_lim$y[2]) {
          main_p_lim$y[2] <- p1_lim$y[2]
        }
        if (p2_lim$y[2] > main_p_lim$x[2]) {
          main_p_lim$x[2] <- p2_lim$y[2]
        }
        p <- p + lims(x = main_p_lim$x, y = main_p_lim$y)
      }
    }
  } else {
    p1 <- ggplot(plotdat, aes(y = x2, fill = level, col = level)) +
      geom_density(alpha = 0.5) +
      common_theme
    p2 <- ggplot(plotdat, aes(x = x1, fill = level, col = level)) +
      geom_density(alpha = 0.5) +
      common_theme
  }
  p1 <- p1 +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal)
  p2 <- p2 +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal)

  p %>%
    aplot::insert_right(p1, width = 0.2) %>%
    aplot::insert_top(p2, height = 0.2) -> p
  p
}


#' 3D plot for b_res
#'
#' @param b_res a b_res object
#' @param Group group vector for color
#' @param metadata metadata contain Group
#' @param ... add
#'
#' @return plotly list
#' @export
#'
#' @examples
#' if (requireNamespace("plotly")) {
#'   data(otutab, package = "pcutils")
#'   b_analyse(otutab, method = "pca", ndim = 3) -> b_res
#'   b_res_3d(b_res, "Group", metadata)
#' }
b_res_3d <- function(b_res, Group, metadata = NULL, ...) {
  lib_ps("plotly", library = FALSE)
  plist <- list()
  if (!is.null(metadata)) {
    if (!Group %in% colnames(metadata)) stop("Group should be one of colnames(metadata)")
    idx <- rownames(metadata) %in% rownames(b_res$sample_site)
    if (!all(rownames(metadata) %in% rownames(b_res$sample_site))) message("rownames don't match in b_res and metadata")
    metadata <- metadata[idx, , drop = FALSE]
    b_res$sample_site <- b_res$sample_site[rownames(metadata), , drop = FALSE]
  } else {
    metadata <- data.frame(row.names = rownames(b_res$sample_site), group = Group)
    Group <- "group"
  }
  for (i in colnames(b_res$sample_eig)[-1]) {
    b_res$sample_site %>% dplyr::select(dplyr::starts_with(i)) -> tmp
    if (ncol(tmp) != 3) next
    plotdat <- data.frame(dim1 = tmp[, 1], dim2 = tmp[, 2], dim3 = tmp[, 3], level = factor(metadata[, Group, drop = TRUE]))
    plotly::plot_ly(plotdat, x = ~dim1, y = ~dim2, z = ~dim3, color = ~level, type = "scatter3d", mode = "markers", ...) %>%
      plotly::layout(title = i) -> plist[[i]]
  }
  return(plist)
}

#' Procrustes Rotation of Two Configurations and PROTEST
#'
#' @param b_res1 Target matrix
#' @param b_res2 Matrix to be rotated
#' @param nperm numbers of permutations to perform
#' @param ... additional
#'
#' @return pro_res
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' b_analyse(otutab, method = "pca") -> b_res1
#' b_analyse(otutab * abs(rnorm(10)), method = "pca") -> b_res2
#' pro_res <- procrustes_analyse(b_res1, b_res2)
#' plot(pro_res, "Group", metadata)
procrustes_analyse <- function(b_res1, b_res2, nperm = 999, ...) {
  if (inherits(b_res1, "b_res")) {
    b_res1 <- b_res1$sample_site[, 2:3]
  }
  message("b_res1 use the ", colnames(b_res1)[1], "\n")
  if (inherits(b_res2, "b_res")) {
    b_res2 <- b_res2$sample_site[, 2:3]
  }
  message("b_res2 use the ", colnames(b_res2)[1], "\n")
  if (nrow(b_res1) != nrow(b_res2)) stop("the row number should be identical for b_res1 and b_res2")
  pro_res <- vegan::protest(b_res1, b_res2, permutations = nperm, ...)
  class(pro_res) <- c("pro_res", class(pro_res))
  pro_res
}

#' Plot pro_res
#'
#' @param x pro_res
#' @param group group
#' @param metadata metadata
#' @param pal pal
#' @param ... add
#'
#' @return a ggplot
#' @exportS3Method
#' @method plot pro_res
#'
plot.pro_res <- function(x, group, metadata = NULL, pal = NULL, ...) {
  pro_res <- x
  X1 <- X2 <- type <- X3 <- X4 <- NULL
  # 提取 Procrustes 分析的坐标
  tab <- cbind(data.frame(pro_res$Yrot), data.frame(pro_res$X))
  colnames(tab) <- c("X1", "X2", "X3", "X4")
  X <- data.frame(pro_res$rotation)

  if (is.null(metadata) && !is.null(group)) {
    md <- data.frame(tab, group = group, check.names = FALSE)
    g_name <- "group"
  } else if (!is.null(metadata) && !is.null(group)) {
    if (!all(rownames(metadata) %in% rownames(tab))) message("rownames dont match in tab and metadata")
    idx <- rownames(metadata) %in% rownames(tab)
    metadata <- metadata[idx, , drop = FALSE]
    tab <- tab[rownames(metadata), , drop = FALSE]
    md <- data.frame(tab, group = metadata[, group, drop = TRUE], check.names = FALSE)
    g_name <- group
  }

  if (is.null(pal)) {
    if (is.continuous(md[, "group"])) {
      pal <- RColorBrewer::brewer.pal(8, "Reds")
    } else {
      pal <- pcutils::get_cols(n = length(unique(md[, "group"])), pal = RColorBrewer::brewer.pal(5, "Set2"))
    }
  }
  point_df <- rbind(
    data.frame(X1 = md$X1, X2 = md$X2, type = "community1", group = md$group),
    data.frame(X1 = md$X3, X2 = md$X4, type = "community2", group = md$group)
  )

  p <- ggplot() +
    geom_point(data = point_df, aes(X1, X2, color = group, shape = type), size = 2) +
    scale_color_manual(values = pal, name = g_name) +
    geom_segment(
      data = md, aes(x = X1, y = X2, xend = X3, yend = X4),
      arrow = arrow(length = unit(0.1, "cm")),
      color = "blue", size = 0.3
    ) +
    theme_classic() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "") +
    geom_vline(xintercept = 0, color = "gray", linetype = 2, size = 0.3) +
    geom_hline(yintercept = 0, color = "gray", linetype = 2, size = 0.3) +
    geom_abline(intercept = 0, slope = X[1, 2] / X[1, 1], size = 0.3) +
    geom_abline(intercept = 0, slope = X[2, 2] / X[2, 1], size = 0.3) +
    labs(subtitle = (paste(expression("M2 ="), round(pro_res$ss, 3), "\n", "P =", pro_res$signif)))
  p
}

match_df <- function(otutab, metadata) {
  if (!setequal(rownames(metadata), colnames(otutab))) message("rownames don't match in tab and metadata")
  idx <- rownames(metadata) %in% colnames(otutab)
  metadata <- metadata[idx, , drop = FALSE]
  otutab <- otutab[, rownames(metadata), drop = FALSE]
  return(list(otutab = otutab, metadata = metadata))
}

#' Permanova between a otutab and a variable
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param envs factors need to test
#' @param norm should normalize?(default:TRUE)
#' @param method adonis/mrpp/anosim/mantel
#' @param dist if use pcoa or nmds, your can choose a dist method (default: bray)
#' @param each test factor one by one, rather than whole
#' @param nperm numbers of permutations to perform
#' @param ... additional
#'
#' @return a g_test object with these columns
#' \item{group}{the test group or factor}
#' \item{r}{relationship}
#' \item{r2}{model R-square}
#' \item{p_value}{model test p_value}
#' \item{sig}{whether significant}
#' @export
#' @references
#' <https://blog.csdn.net/qq_42458954/article/details/110390488>
#' @examples
#' data(otutab, package = "pcutils")
#' permanova(otutab, metadata[, c(2:10)]) -> adonis_res
#' print(adonis_res)
#' plot(adonis_res)
permanova <- function(otutab, envs, norm = TRUE, each = TRUE, method = "adonis",
                      dist = "bray", nperm = 999, ...) {
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

#' Plot g_test
#'
#' @param x a g_test object
#' @param ... add
#'
#' @return ggplot
#' @exportS3Method
#' @method plot g_test
#'
#' @seealso  \code{\link{permanova}}
plot.g_test <- function(x, ...) {
  aa <- x
  group <- r2 <- sig_l <- r <- NULL
  if ("r2" %in% colnames(aa)) {
    aa$group <- factor(aa$group, levels = aa$group[order(aa$r2)]) # 按R2值排序
    aa$sig_l <- cut(aa$p_value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("**", "*", ""))
    p <- ggplot(aa, aes(x = group, y = r2), size = 2) +
      geom_bar(stat = "identity", width = 0.8, color = "black", fill = "#5AAFD1") +
      scale_fill_manual(guide = "none") +
      geom_text(aes(y = r2 * 0.9, label = sig_l), size = 5) + # 可调星号位置
      xlab(NULL) +
      ylab(expression(r^"2")) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_classic()
  }
  if ("r" %in% colnames(aa)) {
    aa$sig_l <- ifelse(aa$sig, "sig", "nosig")
    p <- ggplot(data = aa, aes(x = group, y = r), size = 2) +
      geom_point(aes(size = r, color = sig_l), alpha = 0.8) +
      scale_color_manual(
        values = c("sig" = "red", "nosig" = "blue"),
        breaks = c("sig", "nosig"), name = "Significance", labels = c("p<0.05", "p>0.05")
      ) +
      scale_size_area(max_size = 8) +
      ggrepel::geom_text_repel(label = aa$r) +
      labs(x = NULL, y = "r") +
      geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
      theme_classic()
  }
  return(p)
}


#' Performs graph-based permutation tests
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param metadata metadata
#' @param group one group name in columns of metadata
#' @param nperm numbers of permutations to perform
#' @param ... additional
#'
#' @return ggplot
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("phyloseqGraphTest") && requireNamespace("phyloseq")) {
#'   data(otutab, package = "pcutils")
#'   grap_p_test(otutab, metadata, "Group")
#' }
#' }
grap_p_test <- function(otutab, metadata, group = "Group", nperm = 999, ...) {
  lib_ps("phyloseq", "phyloseqGraphTest", library = FALSE)

  otumat <- otutab
  # random taxon matrix for now. will update to real one later
  taxmat <- matrix(sample(letters, 7 * nrow(otumat), replace = TRUE), nrow = nrow(otumat), ncol = 7)
  rownames(taxmat) <- rownames(otumat)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  sampledata <- phyloseq::sample_data(metadata)
  ps.expo <- phyloseq::phyloseq(
    phyloseq::otu_table(otutab, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxmat), sampledata
  )

  gt1 <- phyloseqGraphTest::graph_perm_test(ps.expo, group, nperm = nperm, ...)
  message("p-value:", gt1$pval, "\n")
  p1 <- phyloseqGraphTest::plot_test_network(gt1)
  p2 <- phyloseqGraphTest::plot_permutations(gt1) +
    labs(subtitle = paste0("p-value = ", gt1$pval, ",   ", nperm, " permutations"))

  return(list(p1, p2))
}


#' Group inter-intra density
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param group group vector
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' gp_dis_density(otutab, metadata["Group"])
gp_dis_density <- function(otutab, group) {
  vegan::vegdist(t(otutab), method = "bray") %>%
    as.b_dist(group_df = group) %>%
    plot.b_dist(., mode = 2)
}


#' RDA
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param env environmental factors
#' @param choose_var should choose variables? use forward step
#' @param norm should normalize? (default:TRUE)
#' @param scale should scale species? (default:FALSE)
#' @param nperm number of permutation
#' @param verbose verbose
#' @param method "rda", "cca", "cap", "dbrda"
#' @param dist The name of the dissimilarity (or distance) index for "cap" or "dbrda", for \code{\link[vegan]{vegdist}}
#' @param direction The direction of the stepwise selection, "both", "forward" or "backward", default is "forward"
#'
#' @return rda/cca
#' @export
#' @seealso \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @examples
#' data(otutab, package = "pcutils")
#' env <- metadata[, 6:10]
#' # RDA
#' myRDA(otutab, env) -> phy.rda
#' RDA_plot(phy.rda, "Group", metadata)
myRDA <- function(otutab, env, norm = TRUE, scale = FALSE, choose_var = FALSE, direction = "forward",
                  nperm = 499, verbose = TRUE, method = "rda", dist = "bray") {
  method <- match.arg(method, c("rda", "cca", "cap", "dbrda"))

  match_res <- match_df(otutab, env)
  otutab <- match_res$otutab
  env <- match_res$metadata

  data.frame(t(otutab)) -> dat
  if (norm) dat.h <- vegan::decostand(dat, "hellinger") else dat.h <- dat
  # env<-decostand(env,method = 'log',MARGIN = 2)#把环境因子进行log转化，以减少同一种环境因子之间本身数值大小造成的影响。

  if (verbose) {
    pcutils::dabiao("Check models", print = TRUE)
    DCA1 <- DCA2 <- NULL
    cat("DCA analysis, select the sorting analysis model according to the first value of the Axis lengths row.
- If it is more than 4.0 - CCA (based on unimodal model, canonical correspondence analysis);
- If it is between 3.0-4.0 - both RDA/CCA;
- If it is less than 3.0 - RDA (based on linear model, redundancy analysis)\n")
    print(vegan::decorana(dat.h) -> dca)
    # p=dca$rproj %>% data.frame() %>%
    #     ggplot(., aes(DCA1, DCA2)) +
    #     labs(title = "DCA") +
    #     geom_point(size = 2) + # 可在这里修改点的透明度、大小
    #     geom_vline(xintercept = 0, color = "gray", size = 0.4) +
    #     geom_hline(yintercept = 0, color = "gray", size = 0.4) +
    #     theme_classic()
    # print(p)
  }

  if (method == "rda") {
    phy.rda <- vegan::rda(formula = dat.h ~ ., data = env, scale = scale)
  } else if (method == "cca") {
    phy.rda <- vegan::cca(formula = dat.h ~ ., data = env, scale = scale)
  } else {
    phy.rda <- vegan::capscale(formula = dat.h ~ ., data = env, scale = scale, distance = dist)
  }

  # print(anova(phy.rda, permutations = how(nperm = 999),by = 'terms'))
  # adonis2(dat.h~.,env,permutations = 9999, method="bray",by = 'terms')
  if (verbose) {
    pcutils::dabiao("Initial Model", print = TRUE)
    cat("Initial cca, vif>20 indicates serious collinearity:\n") # vif>20表示共线性严重。
    print(vegan::vif.cca(phy.rda))
    cat("Initial Model R-square:", (R2a.all <- vegan::RsquareAdj(phy.rda)$adj.r.squared), "\n")
  }

  # 变量的选择
  if (choose_var) {
    # Compare the variance inflation factors
    spe.rda.all <- phy.rda
    # Forward selection using forward.sel()
    # 或者使用ordistep
    if (method == "rda") {
      mod0 <- vegan::rda(formula = dat.h ~ 1, data = env)
    } else if (method == "cca") {
      mod0 <- vegan::cca(formula = dat.h ~ 1, data = env)
    } else {
      mod0 <- vegan::capscale(formula = dat.h ~ 1, data = env, distance = dist)
    }

    pcutils::dabiao("Selecting variables", print = TRUE)
    if (direction == "backward") {
      step.forward <-
        vegan::ordistep(spe.rda.all,
          scope = stats::formula(mod0),
          direction = "backward",
          permutations = permute::how(nperm = nperm)
        )
    } else {
      step.forward <-
        vegan::ordistep(mod0,
          scope = stats::formula(spe.rda.all),
          direction = "forward",
          permutations = permute::how(nperm = nperm)
        )
    }
    if (verbose) {
      ## Parsimonious RDA，简化后的RDA
      pcutils::dabiao("Selected Model", print = TRUE)
      if (is.null(step.forward$CCA)) {
        warning("No variables selected, return the initial model!")
        message("Try to use direction='both', 'backward' or 'forward'")
        return(spe.rda.all)
      }
      # anova(step.forward, permutations = how(nperm = 999),by ='terms' )
      # adonis2(step.forward$call$formula,data = env,permutations = 999, method="bray")
      cat("Selected Model cca, vif>20 means serious collinearity:\n") # vif>20表示共线性严重。
      print(vegan::vif.cca(step.forward))
      cat("Selected Model R-square:", (R2a.pars <- vegan::RsquareAdj(step.forward)$adj.r.squared), "\n")
    }
    phy.rda <- step.forward
  }
  if (verbose) {
    B.sum <- summary(phy.rda)
    pcutils::dabiao("Statistics", print = TRUE)
    cat(B.sum$constr.chi / B.sum$tot.chi, "constrained indicates the degree to which environmental factors explain differences in community structure\n")
    cat(B.sum$unconst.chi / B.sum$tot.chi, "unconstrained means that the environmental factors cannot explain the part of the community structure\n")
  }
  return(phy.rda)
}


#' @export
#' @rdname myRDA
myCCA <- function(otutab, env, norm = TRUE, scale = FALSE, choose_var = FALSE, nperm = 499,
                  verbose = TRUE) {
  myRDA(otutab, env, norm, scale, choose_var, nperm,
    verbose,
    method = "cca"
  )
}

#' @export
#' @rdname myRDA
myCAP <- function(otutab, env, norm = TRUE, scale = FALSE, choose_var = FALSE, nperm = 499,
                  verbose = TRUE, dist = "bray") {
  myRDA(otutab, env, norm, scale, choose_var, nperm,
    verbose,
    method = "cap", dist = dist
  )
}

#' Plot RDA res
#'
#' @param phy.rda rda/cca object
#' @param Group group vector for color
#' @param metadata metadata contain Group
#' @param Group2 mapping point shape
#' @param mode plot mode:1~3
#' @param tri plot variables segments?
#' @param Topn how many variables to show?
#' @param rate segments length rate
#' @param margin plot the margin boxplot?
#' @param box margin plot box or density?
#' @param pal colors for group
#' @param sample_label plot the labels of samples?
#' @param coord_fix fix the coordinates y/x ratio
#' @param bi_text_size biplot text size
#' @param ... add
#' @param stat_ellipse plot the stat_ellipse?
#' @param env_text_param parameters pass to \code{\link[ggplot2]{geom_text}}
#' @param env_rate default 1
#'
#' @seealso  \code{\link{myRDA}}
#' @return ggplot
#' @export
RDA_plot <- function(phy.rda, Group, metadata = NULL, Group2 = NULL, env_rate = 1,
                     mode = 1, tri = FALSE, Topn = 10, rate = 1, margin = FALSE,
                     box = TRUE, pal = NULL, sample_label = TRUE, stat_ellipse = TRUE, coord_fix = FALSE,
                     bi_text_size = 3, env_text_param = NULL, ...) {
  RDA1 <- RDA2 <- contri <- level <- x2 <- x1 <- NULL
  getplotdat <- \(phy.rda, scale = 1){
    if (is.null(phy.rda$CCA)) {
      stop("can be used only with constrained ordination")
    }
    # 提取样方和环境因子排序坐标，前两轴，I 型标尺
    rda.scaling1 <- vegan::scores(phy.rda,
      scaling = scale,
      display = c("sites", "species", "bp")
    )
    rda.site <- data.frame(rda.scaling1$sites)[1:2]
    rda.site$sample <- rownames(rda.site)
    rda.env <- data.frame(rda.scaling1$biplot)[1:2]
    rda.env$sample <- rownames(rda.env)
    rda.spe <- data.frame(rda.scaling1$species)[1:2]
    rda.spe$sample <- rownames(rda.spe)
    return(list(sample_site = rda.site, env_site = rda.env, species_site = rda.spe))
  }

  plotdat <- getplotdat(phy.rda)

  # prepare metadata
  if (!is.null(metadata)) {
    if (!Group %in% colnames(metadata)) stop("Group should be one of colnames(metadata)")
    if (!is.null(Group2)) {
      if (!(Group2 %in% colnames(metadata))) stop("Group2 should be one of colnames(metadata)")
    }
    idx <- rownames(metadata) %in% rownames(plotdat$sample_site)
    if (!all(rownames(metadata) %in% rownames(plotdat$sample_site))) message("rownames don't match in plotdat and metadata")
    metadata <- metadata[idx, , drop = FALSE]
    plotdat$sample_site <- plotdat$sample_site[rownames(metadata), , drop = FALSE]
  } else {
    metadata <- data.frame(row.names = rownames(plotdat$sample_site), group = Group)
    Group <- "group"
    # input group2?
    if (is.null(Group2)) {
      Group2 <- "_shape"
      metadata <- data.frame(metadata, `_group2` = Group2, check.names = FALSE)
      Group2 <- "_group2"
    } else {
      metadata <- data.frame(metadata, group2 = Group2)
      Group2 <- "group2"
    }
  }
  # input group2?
  if (is.null(Group2)) {
    Group2 <- "_shape"
    metadata <- data.frame(metadata, `_group2` = Group2, check.names = FALSE, drop = FALSE)
    Group2 <- "_group2"
  }

  if (is.null(pal)) {
    if (is.numeric(metadata[, Group])) {
      pal <- RColorBrewer::brewer.pal(8, "Reds")
    } else {
      pal <- pcutils::get_cols(n = length(unique(metadata[, Group])), pal = RColorBrewer::brewer.pal(5, "Set2"))
    }
  }

  rda_eig <- (phy.rda$CCA$eig / sum(phy.rda$CCA$eig))[1:2]
  names(rda_eig) -> approach

  colnames(plotdat[[2]])[1:2] <- c("RDA1", "RDA2")
  colnames(plotdat[[3]])[1:2] <- c("RDA1", "RDA2")

  plotdat1 <- data.frame(plotdat[[1]], level = metadata[, Group], level2 = metadata[, Group2])
  colnames(plotdat1)[1:2] <- c("x1", "x2")

  p <- plot_b_like(plotdat1, mode = mode, pal = pal, sample_label = sample_label, stat_ellipse = stat_ellipse, groupname = Group, groupname2 = Group2)

  # envs
  p <- p + labs(
    x = paste(approach[1], ": ", round(100 * rda_eig[1], 2), "%"),
    y = paste(approach[2], ": ", round(100 * rda_eig[2], 2), "%")
  ) +
    geom_segment(
      data = plotdat[[2]], mapping = aes(x = 0, y = 0, xend = RDA1 * env_rate, yend = RDA2 * env_rate),
      arrow = arrow(length = unit(0.2, "cm")), size = 0.5, color = "blue"
    ) +
    do.call(geom_text, pcutils::update_param(list(
      data = plotdat[[2]], mapping = aes(RDA1 * env_rate * 1.1, RDA2 * env_rate * 1.1, label = sample),
      color = "blue", size = 4
    ), env_text_param))
  # geom_text(data = plotdat[[2]], aes(RDA1 * 1.1, RDA2 * 1.1, label = sample),
  #           color = 'blue', size = 4)

  if (tri) {
    plotdat[[3]] %>%
      dplyr::transmute(contri = ((RDA1^2 / rda_eig[1]) + (RDA2^2 / rda_eig[2])) * 100 / 2) %>%
      cbind(plotdat[[3]], .) -> var_contri
    var_contri[, 1:2] * rate -> var_contri[, 1:2]
    var_contri %>% dplyr::top_n(Topn, contri) -> var_contri
    p <- p + ggnewscale::new_scale_color() +
      geom_segment(
        data = var_contri, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = contri),
        size = 0.3, arrow = arrow(length = unit(0.1, "inches"))
      ) +
      geom_text(data = var_contri, aes(
        x = RDA1 * 1.05, y = RDA2 * 1.05,
        color = contri, label = sample
      ), size = 2) +
      scale_color_gradientn(name = "contribution", colours = c("#00AFBB", "#E7B800", "#FC4E07"))
  }

  if (margin) {
    lib_ps("aplot", library = FALSE)
    if (box) {
      p1 <- ggplot(plotdat1, aes(x = level, y = x2, fill = level)) +
        geom_boxplot(outlier.shape = NA) +
        ylab(label = NULL) +
        xlab(label = NULL) +
        scale_fill_manual(values = pal) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = "none",
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = rel(0.8))
        )

      p2 <- ggplot(plotdat1, aes(x = level, y = x1, fill = level)) +
        geom_boxplot(outlier.shape = NA) +
        ylab(label = NULL) +
        xlab(label = NULL) +
        scale_fill_manual(values = pal) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = "none",
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text(size = rel(0.8))
        ) +
        coord_flip()
    } else {
      p1 <- ggplot(plotdat1, aes(y = x2, fill = level, col = level)) +
        geom_density(alpha = 0.5) +
        ylab(label = NULL) +
        xlab(label = NULL) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_classic() +
        theme(
          legend.position = "none",
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()
        )
      p2 <- ggplot(plotdat1, aes(x = x1, fill = level, col = level)) +
        geom_density(alpha = 0.5) +
        ylab(label = NULL) +
        xlab(label = NULL) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_classic() +
        theme(
          legend.position = "none",
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()
        )
    }
    p <- p %>%
      aplot::insert_right(p1, width = 0.2) %>%
      aplot::insert_top(p2, height = 0.2)
  }

  return(p)
}


#' Envfit test for RDA result
#'
#' @param phy.rda a rda result
#' @param env environmental factors
#' @param ... add
#'
#' @return g_test object
#' @export
#' @seealso \code{\link[vegan]{envfit}}
#' @examples
#' data(otutab, package = "pcutils")
#' env <- metadata[, 6:10]
#' # RDA
#' myRDA(otutab, env) -> phy.rda
#' envfitt(phy.rda, env) -> envfit_res
#' plot(envfit_res)
envfitt <- function(phy.rda, env, ...) {
  B.ef <- vegan::envfit(phy.rda, env, ...) # 是做每一个环境因子与群落结构差异的相关性(解释量)
  cor_com <- data.frame(
    group = names(B.ef$vectors$r),
    r2 = B.ef$vectors$r, p_value = B.ef$vectors$pvals
  )
  cor_com$sig <- cor_com$p_value < 0.05
  attributes(cor_com)$method <- "envfit"
  class(cor_com) <- c("g_test", "data.frame")
  return(cor_com)
}
