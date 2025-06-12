# a_diversity==========

#' Calculate a_diversity of otutab
#'
#' @param ... add
#' @param otutab otutab
#'
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' a_diversity(otutab) -> a_res
#' plot(a_res, "Group", metadata)
a_diversity <- function(otutab, ...) {
  UseMethod("a_diversity")
}

#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param method one of "all","richness","chao1","ace","gc","shannon","simpson","pd","pielou","abundance"
#' @param tree a iphylo object match the rownames of otutab
#' @param digits maintance how many digits
#' @rdname a_diversity
#' @return a a_res object
#' @exportS3Method
#' @method a_diversity data.frame
a_diversity.data.frame <- function(otutab, method = c("richness", "shannon"), tree = NULL, digits = 4, ...) {
  all <- c("all", "richness", "chao1", "ace", "gc", "shannon", "simpson", "pd", "pielou", "abundance")
  if (!all(method %in% all)) stop(paste0("methods should be some of ", paste0(all, collapse = ",")))
  if ("all" %in% method) method <- all[-1]
  x <- t(otutab)
  a_res <- data.frame(row.names = colnames(otutab))
  if ("richness" %in% method) {
    Richness <- rowSums(x > 0)
    a_res <- cbind(a_res, Richness)
  }
  if ("abundance" %in% method) {
    Abundance <- rowSums(x)
    a_res <- cbind(a_res, Abundance)
  }
  if ("chao1" %in% method) {
    Chao1 <- vegan::estimateR(x)[2, ]
    a_res <- cbind(a_res, Chao1)
  }
  if ("ace" %in% method) {
    ACE <- vegan::estimateR(x)[4, ]
    a_res <- cbind(a_res, ACE)
  }
  if ("gc" %in% method) {
    Goods_Coverage <- 1 - rowSums(x <= 1) / rowSums(x)
    a_res <- cbind(a_res, Goods_Coverage)
  }
  if ("shannon" %in% method) {
    Shannon <- vegan::diversity(x, index = "shannon", ...)
    a_res <- cbind(a_res, Shannon)
  }
  # 注意，这里是Gini-Simpson 指数
  if ("simpson" %in% method) {
    Simpson <- vegan::diversity(x, index = "simpson", ...)
    a_res <- cbind(a_res, Simpson)
  }
  if ("pielou" %in% method) {
    Pielou_evenness <- vegan::diversity(x, index = "shannon") / log(rowSums(x > 0))
    a_res <- cbind(a_res, Pielou_evenness)
  }
  if ("pd" %in% method) {
    if (is.null(tree)) {
      warning("pd need tree!")
    } else {
      lib_ps("picante", library = FALSE)
      picante::match.phylo.comm(tree, x) -> match_p
      pds <- picante::pd(match_p$comm, match_p$phy, include.root = FALSE)
      PD <- pds[, 1]
      a_res <- cbind(a_res, PD)
      # 净相关指数
      # NRI=-ses.mpd(x,cophenetic(spe_nwk),null.model="taxa.labels")[6]
      # names(NRI) <- 'NRI'
      # 最近邻体指数
      # NTI=-ses.mntd(x,cophenetic(spe_nwk),null.model="taxa.labels")[6]
      # names(NTI) <- 'NTI'
      # result <- cbind(result, PD_whole_tree,NRI,NTI)
    }
  }
  a_res <- round(a_res, digits)
  class(a_res) <- c("a_res", "data.frame")
  return(a_res)
}

#' @param otutab a pc_otu
#'
#' @param method one of "all","richness","chao1","ace","gc","shannon","simpson","pd","pielou"
#' @param tbl which table
#'
#' @exportS3Method
#'
#' @rdname a_diversity
#' @method a_diversity pc_otu
a_diversity.pc_otu <- function(otutab, method = "all", tbl = "otutab", ...) {
  pc <- otutab
  pc_valid(pc)
  otutab <- pc$tbls[[tbl]]
  pc$metas$a_res <- a_diversity.data.frame(otutab, method = method)
  return(pc)
}

#' @param otutab numeric
#' @param ... pass to `a_diversity.data.frame`
#'
#' @exportS3Method
#'
#' @rdname a_diversity
#' @method a_diversity numeric
a_diversity.numeric <- function(otutab, ...) {
  x <- otutab
  return(a_diversity.data.frame(data.frame(Sample = x), ...))
}

#' Plot a_res object
#'
#' @param x a a_res object
#' @param metadata metadata
#' @param group one of colname of metadata
#' @param ... addditional parameters for \code{\link[pcutils]{group_box}} or \code{\link[pcutils]{my_lm}}
#'
#' @return patchwork object,you can change theme with &
#' @exportS3Method
#' @method plot a_res
#'
#' @seealso \code{\link{a_diversity}}
#'
plot.a_res <- function(x, group, metadata, ...) {
  a_res <- x
  a_res <- a_res[rownames(metadata), , drop = FALSE]
  group1 <- metadata[, group]
  if (is.numeric(group1) & !is.factor(group1)) {
    p <- pcutils::my_lm(a_res, group, metadata, ...)
  } else {
    p <- pcutils::group_box(a_res, group, metadata, ...)
  }
  return(p)
}


# test phylogenetic diversity
# if(FALSE){
#   lib_ps("picante")
#   data("phylocom")
#   View(phylocom$sample)
#   ggtree(phylocom$phylo)+geom_tiplab()+theme_tree2()
#   samp=data.frame(a=1:3,b=2:4,c=c(1,2,0),d=c(0,0,3),e=0,row.names = paste0("plot",1:3))
#   read.tree(text = "(c:2,(a:1,b:1):2,d:1,f:1):1;",)->test
#   ggtree(test)+geom_tiplab()+theme_tree2()
#   #prune.sample(samp,test)->test_prune
#   match.phylo.comm(test,samp)->match_p
#   match_p$phy->test;match_p$comm->samp
#
#   test_prune%>%ggtree(.)+geom_tiplab()+theme_tree2()
#   pd(samp,test,include.root = FALSE)
#   pd(samp,test_prune,include.root = FALSE)
#
#   #picante::ses.pd(samp,test)
#   #tree noedes distance
#   cophenetic(test)
#   #每个样方有的物种对的平均谱系距离mpd
#   mpd(samp,cophenetic(test))
#   #随机化mpd
#   mpd(samp,taxaShuffle(cophenetic(test)))
#   #ses.mpd直接计算了
#   ses.mpd(samp,cophenetic(test))->mpd_ses
#   #净相关指数nri,>0聚集，<0发散
#   mpd_ses%>%mutate(nri=-1*(mpd.obs-mpd.rand.mean)/mpd.rand.sd)%>%select(nri)
#   #nti类似nri
#   #mnpd最近谱系距离均值
#   mntd(samp,cophenetic(test))
#   ses.mntd(samp,cophenetic(test))->mnpd_ses
#   mnpd_ses%>%mutate(nti=-1*(mntd.obs-mntd.rand.mean)/mntd.rand.sd)%>%select(nti)
#
#   #beta-mpd
#   comdist(samp,cophenetic(test))
#   #beta-mntd
#   comdistnt(samp,cophenetic(test))
#   #pcd
#   picante::pcd(samp,test)
#   #phylosor
#   picante::phylosor(samp,test)
#   picante::psd(samp,test)
#   picante::raoD(samp,test)
#   picante::unifrac(samp,test)
#
# }


# z_diversity==========
# https://cloud.tencent.com/developer/article/1672945

#' Calculate Zeta Diversity with Distance
#'
#' This function calculates Zeta diversity for each group in the provided otutab.
#'
#' @param otutab A matrix or data frame containing OTU (Operational Taxonomic Unit) counts.
#' @param group_df A data frame containing group information.
#' @param zetadiv_params Additional parameters to be passed to the Zeta.ddecay function from the zetadiv package.
#' @param xy_df Site coordinates.
#'
#' @return zeta_decay
#' @export
#'
#' @examples
#' if (requireNamespace("zetadiv")) {
#'   data(otutab, package = "pcutils")
#'   zeta_decay_result <- z_diversity_decay(otutab, metadata[, c("lat", "long")],
#'     metadata["Group"],
#'     zetadiv_params = list(sam = 10)
#'   )
#'   plot(zeta_decay_result)
#' }
z_diversity_decay <- function(otutab, xy_df, group_df = NULL, zetadiv_params = list()) {
  lib_ps("zetadiv", library = FALSE)

  if (is.null(group_df)) {
    group_df <- data.frame(row.names = colnames(otutab), Group = rep("all", ncol(otutab)), check.names = FALSE)
  }

  zeta_decay <- list()

  for (i in unique(group_df[, 1, drop = TRUE])) {
    tmp_df <- pcutils::trans(pcutils::t2(otutab[, rownames(group_df)[group_df[, 1, drop = TRUE] == i]]), "pa")
    tmp_xy_df <- xy_df[rownames(group_df)[group_df[, 1, drop = TRUE] == i], ]
    tmp_zeta <- do.call(
      zetadiv::Zeta.ddecay,
      update_param(list(
        xy = tmp_xy_df, data.spec = tmp_df, sam = 100, order = 3,
        method.glm = "glm.fit2", confint.level = 0.95,
        normalize = "Jaccard", plot = FALSE
      ), zetadiv_params)
    )
    zeta_decay[[i]] <- tmp_zeta
  }

  class(zeta_decay) <- "zeta_decay"
  return(zeta_decay)
}

#' Plot Zeta Diversity with Distance Results
#'
#' @param x Zeta diversity results obtained from z_diversity_decay function.
#' @param ribbon Logical, whether to add a ribbon to the plot for standard deviation.
#' @param ... Additional arguments to be passed to ggplot2 functions.
#'
#' @return A ggplot object.
#' @exportS3Method
#' @method plot zeta_decay
#'
#' @rdname z_diversity_decay
plot.zeta_decay <- function(x, ribbon = TRUE, ...) {
  zeta_decay <- x
  plot_df <- data.frame()
  distance <- zeta_val <- Group <- fit <- NULL

  for (i in names(zeta_decay)) {
    zeta.bird2 <- zeta_decay[[i]]
    # Predictions
    preds <- stats::predict(zeta.bird2$reg,
      newdata = data.frame(distance.reg = sort(zeta.bird2$distance)),
      type = "link", se.fit = TRUE
    )
    critval <- 1.96

    # Create a data frame for ggplot
    plot_data <- data.frame(
      Group = i,
      distance = sort(zeta.bird2$distance),
      zeta_val = zeta.bird2$zeta.val,
      fit = preds$fit,
      sd = (critval * preds$se.fit)
    )

    plot_df <- rbind(plot_df, plot_data)
  }

  # Plot with ggplot2
  p <- ggplot(plot_df, aes(x = distance, y = zeta_val, color = Group)) +
    geom_point(shape = 16) +
    geom_line(aes(y = fit)) +
    labs(x = "Distance", y = paste0("Zeta diversity (Order ", zeta_decay[[1]]$order, ")"))
  if (ribbon) {
    p <- p +
      geom_ribbon(aes(ymin = fit - sd, ymax = fit + sd, group = Group),
        color = NA, fill = "grey", alpha = 0.5
      )
  }

  return(p)
}

#' Calculate Zeta Diversity
#'
#' This function calculates Zeta diversity for each group in the provided otutab.
#'
#' @param otutab A matrix or data frame containing OTU (Operational Taxonomic Unit) counts.
#' @param group_df A data frame containing group information.
#' @param zetadiv_params Additional parameters to be passed to the Zeta.decline.mc function from the zetadiv package.
#'
#' @return zeta_res
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("zetadiv")) {
#'   data(otutab, package = "pcutils")
#'   zeta_result <- z_diversity(otutab, metadata["Group"], zetadiv_params = list(sam = 10))
#'   plot(zeta_result, lm_model = "exp", text = TRUE)
#' }
#' }
z_diversity <- function(otutab, group_df = NULL, zetadiv_params = list()) {
  lib_ps("zetadiv", library = FALSE)

  if (is.null(group_df)) {
    group_df <- data.frame(row.names = colnames(otutab), Group = rep("all", ncol(otutab)), check.names = FALSE)
  }

  zeta_res <- list()

  for (i in unique(group_df[, 1, drop = TRUE])) {
    tmp_df <- pcutils::trans(pcutils::t2(otutab[, rownames(group_df)[group_df[, 1, drop = TRUE] == i]]), "pa")
    tmp_zeta <- do.call(
      zetadiv::Zeta.decline.mc,
      update_param(list(
        data.spec = tmp_df, orders = 1:5,
        sam = 100, normalize = "Jaccard", plot = FALSE, silent = TRUE
      ), zetadiv_params)
    )
    zeta_res[[i]] <- tmp_zeta
  }

  class(zeta_res) <- "zeta_res"
  return(zeta_res)
}

#' Plot Zeta Diversity Results
#'
#' This function plots the Zeta diversity results obtained from the z_diversity function.
#'
#' @param x Zeta diversity results obtained from z_diversity function.
#' @param lm_model The linear model to be used for fitting ('exp' or 'pl').
#' @param ribbon Logical, whether to add a ribbon to the plot for standard deviation.
#' @param text Logical, whether to add R-squared and p-value text annotations.
#' @param ... Additional arguments to be passed to ggplot2 functions.
#'
#' @return A ggplot object.
#' @exportS3Method
#' @method plot zeta_res
#'
#' @rdname z_diversity
plot.zeta_res <- function(x, lm_model = c("exp", "pl")[1], ribbon = FALSE, text = TRUE, ...) {
  zeta_res <- x

  plot_df <- data.frame()
  p_df <- data.frame()
  `Zeta order` <- `Zeta diversity` <- Group <- V1 <- V2 <- r2 <- NULL
  for (i in names(zeta_res)) {
    zeta.bird2 <- zeta_res[[i]]
    plot_df <- rbind(plot_df, data.frame(
      "Group" = i, "Zeta order" = zeta.bird2$zeta.order,
      "Zeta diversity" = zeta.bird2$zeta.val,
      "sd" = zeta.bird2$zeta.val.sd,
      check.names = FALSE
    ))

    if (lm_model == "exp") {
      tmp_lm <- (zeta.bird2$zeta.exp)
    } else {
      tmp_lm <- (zeta.bird2$zeta.pl)
    }

    p_df <- rbind(p_df, data.frame(
      "Group" = i,
      r2 = round(summary(tmp_lm)$r.squared, 4),
      p = round(anova(tmp_lm)$`Pr(>F)`[1], 4)
    ))
  }

  p <- ggplot(plot_df, aes(x = `Zeta order`, y = `Zeta diversity`, col = Group)) +
    geom_point() +
    geom_line()

  if (ribbon) {
    p <- p + geom_ribbon(aes(ymin = `Zeta diversity` - sd, ymax = `Zeta diversity` + sd, group = Group),
      color = NA, fill = "grey", alpha = 0.5
    )
  }

  lims <- pcutils::ggplot_lim(p)
  p_coor <- pcutils::generate_labels(names(zeta_res), input = c(0.8 * lims$x[2], lims$y[2]), ncols = 1, y_offset = diff(lims$y) * 0.1) %>% as.data.frame()
  p_df <- cbind(p_df, p_coor)

  if (text) {
    p <- p +
      geom_text(data = p_df, aes(x = V1, y = V2, label = paste0("R2= ", r2, "; p= ", p)), show.legend = FALSE)
  }

  return(p)
}
