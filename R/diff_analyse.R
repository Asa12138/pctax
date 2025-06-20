# different analyse==============

#' Difference analysis
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param ctrl the control group, one level of groups
#' @param method one of "deseq2","edger","limma","t.test","wilcox.test"
#' @param log do log transfer for limma?
#' @param add_mini add_mini when calculate the logFC. e.g (10+0.1)/(0+0.1), default `0.5*min(abundance)`
#' @param ... other parameters
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("limma")) {
#'   data(otutab, package = "pcutils")
#'   diff_da(otutab, metadata["Group"], method = "limma") -> res
#'   volcano_p(res)
#'   volcano_p(res, mode = 2)
#' }
#' }
diff_da <- function(otutab, group_df, ctrl = NULL, method = "deseq2", log = TRUE, add_mini = NULL, ...) {
  Group <- group1 <- group2 <- sig <- tax <- NULL
  if (length(method) > 1) {
    all_res <- list()
    for (i in method) {
      pcutils::dabiao("Try ", i)
      res <- tryCatch(
        {
          diff_da(otutab, group_df, ctrl, method = i, log, add_mini)
        },
        error = function(e) {
          NULL
        }
      )
      all_res[[i]] <- res
    }
    return(all_res)
  }
  method <- match.arg(method, c("deseq2", "edger", "limma", "t.test", "wilcox.test"))
  idx <- rownames(group_df) %in% colnames(otutab)
  group_df <- group_df[idx, , drop = FALSE]
  otutab <- otutab[, rownames(group_df), drop = FALSE]
  group_df %>% dplyr::rename(Group = 1) -> meta
  meta$Group <- factor(meta$Group)
  if (!is.null(ctrl)) {
    meta$Group <- relevel(meta$Group, ctrl)
  } else {
    ctrl <- levels(meta$Group)[1]
  }

  if (method == "deseq2") {
    lib_ps("DESeq2", library = FALSE)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = otutab, colData = meta, design = ~Group) # 构建 DESeqDataSet 对象
    dds <- DESeq2::DESeq(dds) # 差异分析

    # results(dds,name = resultsNames(dds)[2])%>%as.data.frame()
    res <- data.frame()
    for (i in 2:length(DESeq2::resultsNames(dds))) {
      DESeq2::results(dds, name = DESeq2::resultsNames(dds)[i]) %>% as.data.frame() -> tmp
      tibble::rownames_to_column(tmp, "tax") -> tmp
      res <- rbind(res, data.frame(tmp, compare = DESeq2::resultsNames(dds)[i]))
    }
    res$compare <- gsub("Group_", "", res$compare) %>% gsub("_vs_", " vs. ", .)
    res$method <- "DESeq2"
  }

  if (method == "edger") {
    lib_ps("edgeR", library = FALSE)
    res <- data.frame()

    for (i in 2:nlevels(meta$Group)) {
      rbind(
        filter(meta, Group == levels(meta$Group)[i]),
        filter(meta, Group == levels(meta$Group)[1])
      ) -> meta1
      group <- meta1$Group
      otutab_f <- otutab[, rownames(meta1)]

      # 数据预处理
      # （1）构建 DGEList 对象
      dgelist <- edgeR::DGEList(counts = otutab_f, group = group)
      # （2）过滤 low count 数据，例如 CPM 标准化（推荐）
      # keep <- rowSums(cpm(dgelist) > 1 ) >= 2
      # dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
      # （3）标准化，以 TMM 标准化为例
      dgelist_norm <- edgeR::calcNormFactors(dgelist, method = "TMM")
      design <- stats::model.matrix(~ factor(meta1$Group))

      # （1）估算基因表达值的离散度
      dge <- edgeR::estimateDisp(dgelist_norm, design, robust = TRUE)
      # （2）模型拟合，edgeR 提供了多种拟合算法
      # 负二项广义对数线性模型
      fit <- edgeR::glmFit(dge, design, robust = TRUE)
      lrt <- edgeR::topTags(edgeR::glmLRT(fit), n = nrow(dgelist$counts))
      lrt$table -> tmp
      tibble::rownames_to_column(tmp, "tax") -> tmp
      res <- rbind(res, data.frame(tmp, compare = paste(levels(meta$Group)[c(i, 1)], collapse = " vs. ")))
    }
    colnames(res) <- c("tax", "log2FoldChange", "logCPM", "LR", "pvalue", "padj", "compare")
    res$method <- "edgeR"
  }

  if (method == "limma") {
    # limma
    lib_ps("limma", library = FALSE)
    otutab_log <- otutab
    if (log) otutab_log <- log(otutab + 1)
    # log后的数据哦

    # df.matrix <- makeContrasts(KO - WT , levels = list)
    res <- data.frame()
    for (i in 2:nlevels(meta$Group)) {
      rbind(
        filter(meta, Group == levels(meta$Group)[i]),
        filter(meta, Group == levels(meta$Group)[1])
      ) -> meta1
      group <- meta1$Group
      otutab_f <- otutab_log[, rownames(meta1)]

      list <- stats::model.matrix(~ 0 + factor(meta1$Group)) # 把group设置成一个model matrix
      colnames(list) <- c("c", "t")
      df.fit <- limma::lmFit(otutab_f, list) ## 数据与list进行匹配

      df.matrix <- limma::makeContrasts(t - c, levels = list)

      fit <- limma::contrasts.fit(df.fit, df.matrix)
      fit <- limma::eBayes(fit)
      tempOutput <- limma::topTable(fit, n = Inf, coef = 1)
      tibble::rownames_to_column(tempOutput, "tax") -> tempOutput
      res <- rbind(res, data.frame(tempOutput, compare = paste(levels(meta$Group)[c(i, 1)], collapse = " vs. ")))
    }
    colnames(res) <- c("tax", "log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B", "compare")
    res$method <- "limma"
  }

  if (method %in% c("t.test", "wilcox.test")) {
    t(otutab) %>%
      data.frame(., check.names = FALSE) %>%
      dplyr::mutate(Group = meta$Group) -> dat
    reshape2::melt(dat, id.vars = "Group") -> dat
    ggpubr::compare_means(
      formula = value ~ Group, data = dat,
      group.by = "variable", method = method,
      p.adjust.method = "fdr", ref.group = ctrl, ...
    ) -> x.all
    x.all <- dplyr::arrange(x.all, group1, group2)
    x.all %>% rename(tax = "variable", "padj" = "p.adj") -> x.all
    x.all$tax <- as.character(x.all$tax)
    x.all$p.format <- as.numeric(x.all$p.format)
    x.all <- dplyr::filter(x.all, group1 == ctrl)
    x.all$compare <- paste0(x.all$group2, " vs. ", x.all$group1)

    pcutils::hebing(otutab, meta$Group) -> tmp
    if (is.null(add_mini)) add_mini <- min(tmp[tmp > 0]) * 0.5
    logfc <- data.frame()
    for (i in 2:nlevels(meta$Group)) {
      logfc <- rbind(logfc, data.frame(
        tax = rownames(tmp),
        compare = paste0(levels(meta$Group)[i], " vs. ", levels(meta$Group)[1]),
        log2FoldChange = log2((tmp[, levels(meta$Group)[i]] + add_mini) / (tmp[, levels(meta$Group)[1]] + add_mini))
      ))
    }
    res <- dplyr::left_join(x.all, logfc)
  }

  res <- data.frame(res)
  res[which(res$log2FoldChange >= 1 & res$padj < 0.05), "sig"] <- "up"
  res[which(res$log2FoldChange <= -1 & res$padj < 0.05), "sig"] <- "down"
  res[which(is.na(res$sig)), "sig"] <- "none"
  res %>% dplyr::mutate(tax1 = ifelse(sig %in% c("up", "down"), tax, "")) -> res

  return(res = res)
}

#' Volcano plot for difference analysis
#'
#' @param res result of `diff_da` which have colnames: tax, log2FoldChange, padj, compare, sig
#' @param logfc log_fold_change threshold
#' @param adjp adjust_p_value threshold
#' @param mode 1:normal; 2:multi_contrast
#' @param number show the tax number
#' @param text text, TRUE
#' @param repel repel, TRUE
#'
#' @return ggplot
#' @export
#'
#' @seealso \code{\link{diff_da}}
volcano_p <- function(res, logfc = 1, adjp = 0.05, text = TRUE, repel = TRUE, mode = 1, number = FALSE) {
  sig <- tax <- log2FoldChange <- padj <- tax1 <- compare <- value <- x <- y <- NULL
  this_method <- unique(res$method)
  if (length(this_method) != 1) stop("Wrong method column")

  res$sig <- NULL
  res[which(res$log2FoldChange >= logfc & res$padj < adjp), "sig"] <- "up"
  res[which(res$log2FoldChange <= -logfc & res$padj < adjp), "sig"] <- "down"
  res[which(is.na(res$sig)), "sig"] <- "none"
  if (!"tax1" %in% colnames(res)) res %>% dplyr::mutate(tax1 = ifelse(sig %in% c("up", "down"), tax, "")) -> res
  res <- dplyr::arrange(res, -log2FoldChange)
  res$label <- ifelse(res$sig %in% c("up", "down"), "Sig", "Non-sig")
  cols <- setNames(c("#2f5688", "#BBBBBB", "#CC0000"), c("down", "none", "up"))

  if (number) {
    label <- dplyr::count(res, sig) %>% dplyr::mutate(label = paste0(sig, ": ", n))
    res$sig <- setNames(label$label, label$sig)[res$sig]
    cols <- setNames(cols, setNames(label$label, label$sig)[names(cols)])
  }

  if (mode == 1) {
    # unique(res$compare)
    # 一对比较的火山图
    res %>% dplyr::filter(!is.na(padj)) -> dat

    pp <- ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point() +
      scale_color_manual(values = cols, na.value = "#BBBBBB") + # 确定点的颜色
      facet_wrap(. ~ compare, scales = "free") +
      labs(title = this_method) +
      theme_bw(base_size = 14) + # 修改图片背景
      theme(
        legend.title = element_blank() # 不显示图例标题
      ) +
      ylab("-log10 (p-adj)") + # 修改y轴名称
      xlab("log2 (FoldChange)") + # 修改x轴名称
      geom_vline(xintercept = c(-logfc, logfc), lty = 3, col = "black", lwd = 0.5) + # 添加垂直阈值|FoldChange|>2
      geom_hline(yintercept = -log10(adjp), lty = 3, col = "black", lwd = 0.5) # 添加水平阈值padj<0.05

    if (text) {
      if (repel) {
        pp <- pp + ggrepel::geom_text_repel(aes(label = tax1), size = 2, show.legend = FALSE)
      } else {
        pp <- pp + geom_text(aes(label = tax1), size = 2, show.legend = FALSE)
      }
    }
  }
  if (mode == 2) {
    # 多对比较的火山图
    res %>%
      dplyr::filter(abs(log2FoldChange) > logfc) %>%
      dplyr::filter(!is.na(padj)) -> dat
    res %>%
      dplyr::group_by(compare) %>%
      dplyr::summarise(max(log2FoldChange)) %>%
      reshape2::melt() -> bardf
    res %>%
      dplyr::group_by(compare) %>%
      dplyr::summarise(min(log2FoldChange)) %>%
      reshape2::melt() -> bardf1

    p1 <- ggplot() +
      geom_bar(
        data = bardf,
        mapping = aes(x = compare, y = value), stat = "identity",
        fill = "#dcdcdc", alpha = 0.6
      ) +
      geom_bar(
        data = bardf1,
        mapping = aes(x = compare, y = value), stat = "identity",
        fill = "#dcdcdc", alpha = 0.6
      )
    p2 <- p1 + geom_jitter(
      data = dat,
      aes(x = compare, y = log2FoldChange, color = label),
      size = 1.2,
      width = 0.4
    )

    if (text) {
      if (repel) {
        p2 <- p2 + ggrepel::geom_text_repel(data = filter(dat, tax1 != ""), aes(x = compare, y = log2FoldChange, label = tax1), col = "red", size = 3, force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last"))
      } else {
        p2 <- p2 + geom_text(data = filter(dat, tax1 != ""), aes(x = compare, y = log2FoldChange, label = tax1), col = "red", size = 3)
      }
    }

    coldf <- data.frame(x = unique(res$compare), y = 0)
    p3 <- p2 + geom_tile(
      data = coldf,
      aes(x = x, y = y, fill = x),
      height = 0.4,
      color = "black",
      alpha = 0.6,
      show.legend = FALSE
    ) +
      labs(x = "Compares", y = "log2 (FoldChange)") +
      geom_text(
        data = coldf,
        aes(x = x, y = y, label = x),
        size = 6,
        color = "white"
      ) +
      scale_color_manual(
        name = NULL,
        values = c("black", "red")
      )


    pp <- p3 + theme_minimal(base_size = 14) +
      # coord_cartesian(ylim = c(-2,4))+
      theme(
        axis.title = element_text(size = 13, color = "black", face = "bold"),
        axis.line.y = element_line(color = "black", size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        # legend.direction = "h",
        legend.text = element_text(size = 15)
      )
  }
  return(pp)
}

#' Difference analysis
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param mode 1~2
#' @param text_df text_df
#' @param text_x text_x
#' @param text_angle text_angle
#' @param errorbar top, bottom, none
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' multi_bar(otutab[1:10, ], metadata["Group"])
multi_bar <- function(otutab, group_df, mode = 1, text_df = NULL,
                      text_x = NULL, text_angle = -90, errorbar = "bottom") {
  Group <- tax <- abundance <- se <- label <- NULL
  idx <- intersect(rownames(group_df), colnames(otutab))
  group_df <- group_df[idx, , drop = FALSE]
  otutab <- otutab[, idx, drop = FALSE]
  group_df %>% dplyr::rename(Group = 1) -> meta
  # meta$Group=factor(meta$Group)

  t(otutab) %>%
    data.frame(., check.names = FALSE) %>%
    mutate(Group = meta$Group) -> dat
  reshape2::melt(dat, id.vars = "Group", variable.name = "tax", value.name = "abundance") -> dat
  dat$tax <- factor(dat$tax, levels = rownames(otutab))

  low <- min(dat$abundance)
  high <- max(dat$abundance)

  if (mode == 1) {
    meandf <- group_by(dat, Group, tax) %>% summarise(mean = mean(abundance), sd = sd(abundance), se = sd(abundance) / sqrt(n() - 1))
    p <- ggplot(data = meandf, aes(y = tax, x = mean))

    if (errorbar == "bottom") {
      p <- p +
        geom_errorbar(aes(xmin = mean - se, xmax = mean + se, group = Group),
          stat = "identity",
          position = position_dodge(width = 0.9), width = 0.25
        ) +
        geom_bar(aes(fill = Group),
          stat = "identity", position = position_dodge(width = 0.9),
          width = 0.8, size = 0.2, color = "black"
        )
    } else if (errorbar == "top") {
      p <- p +
        geom_bar(aes(fill = Group),
          stat = "identity", position = position_dodge(width = 0.9),
          width = 0.8, size = 0.2, color = "black"
        ) +
        geom_errorbar(aes(xmin = mean - se, xmax = mean + se, group = Group),
          stat = "identity",
          position = position_dodge(width = 0.9), width = 0.25
        )
    } else {
      p <- p +
        geom_bar(aes(fill = Group),
          stat = "identity", position = position_dodge(width = 0.9),
          width = 0.8, size = 0.2, color = "black"
        )
    }

    if (is.null(text_x)) text_x <- -0.05 * (high)
  }
  if (mode == 2) {
    p <- ggplot(data = dat, aes(y = tax, x = abundance)) +
      geom_boxplot(aes(color = Group), outlier.shape = NA)
    if (is.null(text_x)) text_x <- low - 0.05 * (high - low)
  }
  if (!is.null(text_df)) {
    text_df <- text_df[rownames(otutab), , drop = FALSE]
    text_df$tax <- rownames(text_df)
    colnames(text_df)[1] <- "label"
    p <- p + geom_text(data = text_df, mapping = aes(x = text_x, y = tax, label = label), angle = text_angle)
  }
  p + pctax_theme +
    theme(legend.position = "top")
  # lib_ps("ggpubr")
  # lib_ps("ggfun")
  # ggpubr::ggbarplot(dat,y="abundance",x="tax",fill = "Group",add = "mean_se",position = position_dodge(width = 0.8))+coord_flip()+
  #   ggfun::theme_stamp(axis = "x")
}

#' Difference analysis
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' multi_conf(otutab[1:10, 1:12], metadata["Group"])
multi_conf <- function(otutab, group_df) {
  Group <- tax <- abundance <- se <- label <- NULL
  diff_mean <- conf_lower <- conf_upper <- p.value <- NULL

  idx <- intersect(rownames(group_df), colnames(otutab))
  group_df <- group_df[idx, , drop = FALSE]
  otutab <- otutab[, idx, drop = FALSE]
  group_df %>% dplyr::rename(Group = 1) -> meta
  meta$Group <- factor(meta$Group)

  # meta$Group=factor(meta$Group)
  apply(otutab, 1, \(i)t.test(i ~ meta$Group)) -> a
  lapply(a, \(i)c(-diff(i$estimate), i$conf.int[1], i$conf.int[2], i$p.value) %>% unname()) -> conf_res
  conf_res_df <- data.frame(tax = names(conf_res), t2(data.frame(conf_res)))
  colnames(conf_res_df) <- c("tax", "diff_mean", "conf_lower", "conf_upper", "p.value")
  conf_res_df$tax <- factor(conf_res_df$tax, levels = rownames(otutab))
  conf_res_df$Group <- ifelse(conf_res_df$diff_mean > 0, levels(meta$Group)[1], levels(meta$Group)[2])

  p1 <- ggplot(conf_res_df, aes(x = diff_mean, y = tax)) +
    geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper, group = Group), width = 0.2) +
    geom_point(aes(fill = Group), shape = 21, size = 4) +
    geom_vline(xintercept = 0, linetype = 2) +
    pctax_theme +
    theme(legend.position = "top")
  p2 <- ggplot(conf_res_df, aes(x = 0, y = tax, label = format.pval(p.value))) +
    geom_text() +
    theme_void()
  patchwork::wrap_plots(p1, p2, ncol = 2)
}


#' Stamp style plot
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param set_order set order of factor levels
#' @param pal palette
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' if (requireNamespace("ggfun")) stamp_plot(otutab[1:10, 1:12], metadata["Group"])
stamp_plot <- function(otutab, group_df, set_order = NULL, pal = NULL) {
  lib_ps("ggfun", library = FALSE)

  # otutab=trans(otutab,"total")*100
  p1 <- multi_bar(otutab, group_df)
  p1 <- p1 + ggfun::theme_stamp(colour = c("white", "grey90")) +
    theme(panel.grid.major.x = element_blank(), text = element_text(face = "bold")) +
    labs(y = NULL, x = "Mean proportion (%)")

  if (!is.null(set_order)) p1$data$tax <- change_fac_lev(p1$data$tax, levels = set_order)

  if (is.null(pal)) pal <- c("red", "blue")
  if (is.null(names(pal))) names(pal) <- levels(factor(p1$data$Group))

  p <- multi_conf(otutab, group_df)
  if (!is.null(set_order)) {
    p[[1]]$data$tax <- change_fac_lev(p[[1]]$data$tax, levels = set_order)
    p[[2]]$data$tax <- change_fac_lev(p[[2]]$data$tax, levels = set_order)
  }
  p2 <- p[[1]] + ggfun::theme_stamp(colour = c("white", "grey90")) + theme(panel.grid.major.x = element_blank()) +
    labs(y = NULL, x = "Difference in mean proportions (%)", title = "95% confidence intervals") +
    theme(
      text = element_text(face = "bold"),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 15, face = "bold", colour = "black", hjust = 0.5)
    )
  p3 <- p[[2]] + ggtitle("P-values") + theme(
    text = element_text(face = "bold"),
    plot.title = element_text(size = 15, face = "bold", colour = "black", hjust = 0.5)
  )
  p3$layers[[1]]$aes_params$fontface <- "bold"
  comp <- patchwork::wrap_plots(p1, p2, p3) + patchwork::plot_layout(widths = c(4, 6, 2))
  return(comp & scale_fill_manual(values = pal))
}


#' Get mean and type
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @return No value
#' @export
get_diff_type <- function(otutab, group_df) {
  idx <- rownames(group_df) %in% colnames(otutab)
  group_df <- group_df[idx, , drop = FALSE]
  otutab <- otutab[, rownames(group_df), drop = FALSE]

  group <- group_df[[1]] %>% as.factor()
  pcutils::hebing(otutab, group) -> tmp
  apply(tmp, 1, function(a) which(a == max(a))[[1]]) -> tmp$group
  apply(tmp, 1, function(a) {
    return(colnames(tmp)[a[[ncol(tmp)]]])
  }) -> tmp$Type
  tmp$Type <- factor(tmp$Type, levels = levels(group))
  tmp$tax <- rownames(tmp)
  tmp
}


#' KW test
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param method "kruskal.test", see \code{\link[ggpubr]{compare_means}}
#'
#' @return res
#' @export
#'
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' kwtest(otutab, metadata["Group"]) -> res
#' bbtt(res, pvalue = "p.format")
#' }
kwtest <- function(otutab, group_df, method = "kruskal.test") {
  group <- group_df[[1]] %>% as.factor()
  t(otutab) %>%
    data.frame(., check.names = FALSE) %>%
    dplyr::mutate(Group = group) -> dat
  reshape2::melt(dat, id.vars = "Group") -> dat
  ggpubr::compare_means(
    formula = value ~ Group, data = dat, group.by = "variable",
    method = method, p.adjust.method = "fdr"
  ) -> x.all
  x.all %>% rename(tax = "variable") -> x.all
  x.all$tax <- as.character(x.all$tax)
  x.all$p.format <- as.numeric(x.all$p.format)

  tmp <- get_diff_type(otutab, group_df)
  x.all <- dplyr::left_join(x.all, tmp, by = "tax")
  return(x.all)
}

#' ALDEX
#'
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param otutab otutab
#'
#' @return diff
#' @export
#' @references
#' <https://cloud.tencent.com/developer/article/1621879>
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ALDEx2")) {
#'   data(otutab, package = "pcutils")
#'   ALDEX(otutab, metadata["Group"]) -> res
#'   res %>%
#'     dplyr::top_n(9, -glm.eBH) %>%
#'     .[, "tax"] -> sig
#'   data.frame(t(otutab[sig, ])) %>% pcutils::group_box(., "Group", metadata)
#' }
#' }
ALDEX <- function(otutab, group_df) {
  lib_ps("ALDEx2", library = FALSE)
  we.eBH <- glm.eBH <- NULL
  group <- group_df[[1]] %>% as.factor()
  if (nlevels(factor(group)) == 2) {
    x.all <- ALDEx2::aldex(otutab, group,
      mc.samples = 16, test = "t", effect = TRUE,
      include.sample.summary = FALSE, denom = "all", verbose = FALSE, paired.test = FALSE
    )
    data.frame(tax = rownames(x.all), x.all) -> x.all
    x.all %>% mutate(p.signif = ifelse(we.eBH >= 0.05, "ns", ifelse(we.eBH >= 0.01, "*", ifelse(we.eBH > 0.001, "**", "***")))) -> x.all
  }
  if (nlevels(factor(group)) > 2) {
    x.all <- ALDEx2::aldex(otutab, group,
      mc.samples = 16, test = "kw", effect = TRUE,
      include.sample.summary = FALSE, denom = "all", verbose = FALSE, paired.test = FALSE
    )
    data.frame(tax = rownames(x.all), x.all) -> x.all
    x.all %>% mutate(p.signif = ifelse(glm.eBH >= 0.05, "ns", ifelse(glm.eBH >= 0.01, "*", ifelse(glm.eBH > 0.001, "**", "***")))) -> x.all
  }

  tmp <- get_diff_type(otutab, group_df)
  x.all <- dplyr::left_join(x.all, tmp, by = "tax")
  return(x.all)
}

# ANCOM2<-function(otutab,group_df){
#   res = ANCOM(otutab, group_df, struc_zero = NULL,main_var = colnames(group_df))
#   res$out->out
#   rownames(out)<-out$taxa_id
#   hebing(otutab,metadata[,group])->tmp
#   apply(tmp, 1, function(a)which(a==max(a))[[1]])->tmp$group
#   apply(tmp, 1, function(a){return(colnames(tmp)[a[[4]]]) })->tmp$Type
#   rename(out,tax=taxa_id)->out
#   cbind(out,tmp[out$tax,])->out
#   cbind(out,otutab[out$tax,])->out
#
#   # Step 3: Volcano Plot
#   n_taxa = nrow(otutab)
#   # Cutoff values for declaring differentially abundant taxa
#   cut_off = c(0.9 * (n_taxa - 1), 0.8 * (n_taxa - 1), 0.7 * (n_taxa - 1), 0.6 * (n_taxa - 1))
#   names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")
#   # Annotation data
#   dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
#   res$fig = res$fig +
#     geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
#     geom_text(data = dat_ann, aes(x = x, y = y, label = label),
#               size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
#   return(res)
# }


#' ggdotchart for diff analysis
#'
#' @param res result of ALDEX or kwtest
#' @param pvalue the name of pvaule
#' @param topN topN
#'
#' @return ggplot
#' @export
bbtt <- function(res, pvalue = "glm.eBH", topN = 20) {
  Type <- tax <- NULL
  res %>%
    dplyr::arrange(get(pvalue)) %>%
    head(topN) %>%
    dplyr::mutate(group = Type, tax = factor(tax, levels = rev(tax))) -> pres1
  p1 <- ggpubr::ggdotchart(pres1,
    x = "tax", y = pvalue, color = "group",
    dot.size = 7, # palette = classcol, # 修改颜色
    sorting = "ascending",
    add = "segments",
    add.params = list(color = "lightgray", size = 1.5), # 添加棒子
    label = "Type", # 添加label
    font.label = list(color = "white", size = 8, vjust = 0.5),
    rotate = TRUE,
    ggtheme = pctax_theme # 改变主题
  ) + labs(x = "")
  p1
  # pres1[,levels(pres1$Type)]%>%pcutils::trans(.,1,method = 'normalize')%>%
  #   dplyr::mutate(tax=rownames(.),tax=factor(tax,levels = rev(tax)))->tmp
  # p2<-reshape2::melt(tmp)%>%
  #   ggplot(.,aes(x=variable,y=tax,fill=value))+geom_raster()+
  #   scale_fill_gradientn(colours = c("#1789E6", "#FFFFFF", "#B3192B"))+
  #   theme_void()+
  #   theme(axis.text.x = element_text(),legend.position = 'top')
  # lib_ps("patchwork",library = FALSE)
  # p2+p1+plot_layout(widths = c(1.5,3))
}

#' RandomForest
#'
#' @param otutab otutab
#' @param topN default: 10
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return diff
#' @export
#'
#' @examples
#' if (requireNamespace("randomForest")) {
#'   data(otutab, package = "pcutils")
#'   suijisenlin(otutab, metadata["Group"]) -> rf_res
#' }
suijisenlin <- function(otutab, group_df, topN = 10) {
  group <- group_df[[1]]
  MeanDecreaseAccuracy <- tax <- NULL
  lib_ps("randomForest", library = FALSE)
  t(otutab) %>%
    data.frame() %>%
    mutate(Group = group) %>%
    randomForest::randomForest(Group ~ ., data = ., importance = TRUE, proximity = TRUE) -> randomforest

  # plot(randomforest)
  # 查看变量的重要性
  randomforest$importance %>%
    data.frame() %>%
    tibble::rownames_to_column("tax") -> im

  im %>%
    dplyr::top_n(topN, wt = MeanDecreaseAccuracy) %>%
    dplyr::arrange(-MeanDecreaseAccuracy) %>%
    dplyr::mutate(tax = factor(tax, levels = rev(tax))) -> mim

  imp <- ggpubr::ggscatter(mim,
    x = "MeanDecreaseAccuracy", y = "tax", color = "#7175E3",
    size = "MeanDecreaseAccuracy", ggtheme = pctax_theme
  ) + ylab("") +
    theme(legend.position = "none")

  # pcutils::hebing(otutab,group =group)->tmp
  # tmp[as.character(mim$tax),]%>%decostand(.,1,method = 'normalize')%>%
  #   mutate(tax=rownames(.),tax=factor(tax,levels = rev(tax)))->tmp
  #
  # p2<-melt(tmp)%>%
  #   ggplot(.,aes(x=variable,y=tax,fill=value))+geom_raster()+
  #   scale_fill_gradientn(colours = c("#003366","white","#990033"))+theme_pubr(legend = 'right')+
  #   theme_void()+
  #   theme(axis.text.x = element_text(),legend.position = 'right')
  # imp<-imp+p2+plot_layout(widths = c(2.5,1.5))
  return(list(im = im, imp = imp))
}

#' Time series analysis
#'
#' @param otu_time otutab hebing by a time variable
#' @param n_cluster number of clusters
#' @param min.std min.std
#'
#' @return time_cm
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   data(otutab, package = "pcutils")
#'   otu_time <- pcutils::hebing(otutab, metadata$Group)
#'   time_by_cm(otu_time, n_cluster = 4) -> time_cm_res
#'   plot(time_cm_res)
#' }
#' }
time_by_cm <- function(otu_time, n_cluster = 6, min.std = 0) {
  lib_ps("Mfuzz", "methods", library = FALSE)

  otu_time <- as.matrix(otu_time)

  mfuzz_class <- methods::new("ExpressionSet", exprs = otu_time)
  mfuzz_class <- Mfuzz::filter.NA(mfuzz_class, thres = 0.25)
  mfuzz_class <- Mfuzz::fill.NA(mfuzz_class, mode = "mean")
  mfuzz_class <- Mfuzz::filter.std(mfuzz_class, min.std = min.std)
  # 标准化数据
  mfuzz_class <- Mfuzz::standardise(mfuzz_class)
  mfuzz_cluster <- Mfuzz::mfuzz(mfuzz_class, c = n_cluster, m = Mfuzz::mestimate(mfuzz_class))
  # mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = rownames(otu_time))

  cbind(mfuzz_class@assayData$exprs %>% as.data.frame(),
    id = mfuzz_class@assayData$exprs %>% as.data.frame() %>% rownames(),
    cluster = mfuzz_cluster$cluster,
    membership = mfuzz_cluster$membership %>% apply(., 1, max)
  ) -> plotdat
  class(plotdat) <- c("time_cm", "data.frame")
  pcutils::del_ps("Mfuzz")
  return(plotdat)
}

#' Plot time_cm
#'
#' @param x time_cm
#' @param mem_thr membership threshold
#' @param ... add
#'
#' @return ggplot
#' @exportS3Method
#' @method plot time_cm
#' @rdname c_means
plot.time_cm <- function(x, mem_thr = 0.6, ...) {
  membership <- value <- NULL
  fancy.blue <- c(
    c(255:0), rep(0, length(c(255:0))),
    rep(0, length(c(255:150)))
  )
  fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
  fancy.red <- c(
    c(0:255), rep(255, length(c(255:0))),
    c(255:150)
  )
  colo <- grDevices::rgb(b = fancy.blue / 255, g = fancy.green / 255, r = fancy.red / 255)

  plotdat <- x
  plotdat %>% filter(membership >= mem_thr) -> plotdat1
  plotdat <- reshape2::melt(plotdat1, id.vars = c("id", "cluster", "membership"), variable.name = "time")
  plotdat$cluster <- factor(plotdat$cluster)

  ggplot(data = plotdat, aes(x = time, y = value, group = id, col = membership, alpha = membership)) +
    facet_wrap(cluster ~ .) +
    scale_color_gradientn(colours = colo) +
    geom_line(size = 0.8) +
    pctax_theme +
    scale_x_discrete(expand = c(0, 0.2)) +
    theme(plot.margin = unit(c(1, 2, 1, 1), "lines"))
}


# mpse_da<-function(otutab,metadata_g,taxonomy,alpha=0.05){
#   lib_ps("MicrobiotaProcess")
#   data.frame(tax=taxonomy%>%apply(., 1, \(x)paste(unlist(x),collapse = "|")),otutab,check.names = FALSE)->motu
#   write.table(motu,file = "./tmp",row.names = FALSE,sep = "\t",quote = FALSE)
#   if(ncol(metadata_g)!=2)stop("metadata_g need two columns, first is id, second is group")
#   colnames(metadata_g)[2]<-"Group"
#   MicrobiotaProcess::mp_import_metaphlan(profile="./tmp", mapfilename=metadata_g)->mpse
#   file.remove("./tmp")
#
#   mpse%>%
#     mp_diff_analysis(
#       .abundance = Abundance,
#       .group = Group,
#       first.test.alpha = alpha
#     )->mpse2
#
#   mpse2%<>%
#     mp_cal_abundance( # for each samples
#       .abundance = RareAbundance
#     ) %>%
#     mp_cal_abundance( # for each groups
#       .abundance=RareAbundance,
#       .group=Group
#     )
#
#   taxa.tree <- mpse2 %>%
#     mp_extract_tree(type="taxatree")
#   #taxa.tree %>% dplyr::select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Group, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))
#
#   p1<-mpse2%>%
#     mp_plot_diff_res(
#       group.abun = TRUE,
#       pwidth.abun=0.1,
#       tiplab.size=1
#     )
#
#   tibble::as_tibble(taxa.tree)%>%filter(!is.na(LDAmean))->saa1
#
#   tidyr::unnest(saa1,RareAbundanceBySample)->saa2
#
#   pp1<-ggplot(saa2, aes(x=label, y=RelRareAbundanceBySample,fill=Group)) +
#     geom_boxplot(aes(),width = 0.5) + ylab(label = "Abundance")+
#     #scale_fill_d3()+
#     coord_flip()+
#     ggfun::theme_stamp(
#       colour = c('grey90', 'white'),
#       axis = 'x',
#       axis.line.x = element_line(),
#       axis.title.y = element_blank(),legend.direction = "vertical")
#
#   pp2<-ggplot(saa1, aes(label, LDAmean)) +
#     ylab("LDA SCORE (log 10)")  +
#     geom_point(aes(col = Sign_Group))+
#     #scale_color_d3()+
#     #geom_bar(stat = "identity", aes(fill = Sign_Group),width = 0.5) + scale_fill_d3()+
#     coord_flip()+
#     ggfun::theme_stamp(
#       colour = c('grey90', 'white'),
#       axis = 'x',
#       axis.line.x = element_line(),
#       axis.title.y = element_blank(),
#       axis.line.y= element_blank(),
#       axis.ticks.y  = element_blank(),
#       axis.text.y=element_blank(),
#       legend.position = "none")
#   lib_ps("aplot",library = FALSE)
#   p2<-pp1%>%insert_right(pp2,width = 0.6)
#
#   p3<-mpse2%>%
#     mp_plot_diff_cladogram(
#       label.size = 2.5,
#       hilight.alpha = .3,
#       bg.tree.size = .5,
#       bg.point.size = 2,
#       bg.point.stroke = .25
#     )
#
#   detach("package:MicrobiotaProcess")
#   return(list(tree=taxa.tree,p1=p1,p2=p2,p3=p3))
# }
