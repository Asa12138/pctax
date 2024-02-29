# rare==========

#' Rarefy a otutab
#'
#' @param otutab otutab
#' @param sample number
#'
#' @return a rarefied otutab
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' rarefaction(otutab)
rarefaction <- function(otutab, sample = NULL) {
    lib_ps("vegan", library = FALSE)
    if (is.null(sample)) sample <- min(colSums(otutab))
    otutab1 <- vegan::rrarefy(t(otutab), sample = sample)
    return(as.data.frame(t(otutab1)))
}

#' Rare the sample
#'
#' @param otutab otutab
#' @param rep repeats number
#' @param count_cutoff cutoff to be 0
#'
#' @export
#' @return ggplot
#' @examples
#' data(otutab, package = "pcutils")
#' a <- rare_curve_sample(otutab)
#' plot(a)
rare_curve_sample <- function(otutab, rep = 30, count_cutoff = 1) {
    # 默认值绘制样本箱线图稀释曲线
    count <- otutab
    count[count < count_cutoff] <- 0
    count[count >= count_cutoff] <- 1
    n <- dim(count)[2]
    x <- unique(as.integer(seq(1, n, length.out = n)))
    result <- data.frame(sample = c(0), richness = c(0))
    for (i in x) {
        m <- choose(n, i)
        if (m > rep) {
            m <- rep
        }
        for (j in 1:m) {
            idx <- sample(n, i)
            temp <- count[, idx, drop = FALSE]
            temp1 <- temp[rowSums(temp) > 0, , drop = FALSE]
            result <- rbind(result, c(i, dim(temp1)[1]))
        }
    }
    result <- result[-1, ]
    result$sample <- factor(result$sample, levels = unique(result$sample))
    attributes(result)$type <- "sample"
    class(result) <- c("rare_res", class(result))
    return(result)
}

#' Plot a rare curve
#'
#' @param x AOR object
#' @param ... add
#'
#' @rdname rare_curve_sample
#' @method plot rare_res
#' @exportS3Method
#' @return ggplot
plot.rare_res <- function(x, ...) {
    result <- x
    richness <- rare <- NULL
    if (attributes(result)$type == "sample") {
        p <- ggplot(result, aes(x = sample, y = richness, fill = sample)) +
            geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.5, size = 0.1) +
            xlab("Sample size") +
            ylab("Richness") +
            theme(legend.position = "none") +
            theme_bw()
    }
    if (attributes(result)$type == "species") {
        p <- ggplot(result, aes(x = as.numeric(rare), y = alpha, color = sample)) +
            geom_line() +
            geom_errorbar(aes(ymin = alpha - sd, ymax = alpha + sd)) +
            geom_point(size = 1) +
            labs(
                x = ifelse(attributes(result)$mode == 1, "reads (count)", "reads (cpm)"),
                y = attributes(result)$method
            ) +
            theme_bw()
    }
    p
}

#' Rare the species
#'
#' @param otutab otutab
#' @param step default 2000
#' @param method one of "richness","chao1","ace","gc","shannon","simpson","pd","pielou"
#' @param mode 1 for little table, 2 for big
#' @param verbose verbose
#' @param threads use how many threads to calculate (default:1)
#' @param reps reps
#'
#' @export
#' @return ggplot
#' @examples
#' data(otutab, package = "pcutils")
#' a <- rare_curve_species(otutab, mode = 1)
#' plot(a)
rare_curve_species <- function(otutab, step = 2000, method = "richness", mode = 2,
                               reps = 3, threads = 1, verbose = TRUE) {
    i <- rare <- NULL
    # 根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
    alpha_curves <- \(x, step, method = "richness", tree = NULL) {
        x_nrow <- nrow(x)
        rare <- rowSums(x)

        alpha_rare <- data.frame()
        for (i in 1:x_nrow) {
            step_num <- seq(0, rare[i], step)
            if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
            alpha_rare_i <- NULL
            for (step_num_n in step_num) {
                pctax::a_diversity(data.frame(t(vegan::rrarefy(x[i, ], step_num_n))), method = method, tree = tree) -> tmp_i
                alpha_rare_i <- c(alpha_rare_i, unlist(tmp_i))
            }
            alpha_rare <- rbind(
                alpha_rare,
                data.frame(rare = step_num, alpha = alpha_rare_i, sample = rownames(x)[i])
            )
        }
        alpha_rare
    }

    if (mode == 1) {
        t(otutab) -> totu
    } else {
        t(round(pcutils::trans(otutab, "cpm"), 0)) -> totu
    }

    # parallel
    # main function
    loop <- function(i) {
        alpha_curves(totu, step = step, method = "richness")
    }
    {
        if (threads > 1) {
            pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
            if (verbose) {
                pb <- utils::txtProgressBar(max = reps, style = 3)
                opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
            } else {
                opts <- NULL
            }
            cl <- snow::makeCluster(threads)
            doSNOW::registerDoSNOW(cl)
            res <- foreach::`%dopar%`(
                foreach::foreach(i = seq_len(reps), .options.snow = opts),
                loop(i)
            )
            snow::stopCluster(cl)
            gc()
        } else {
            res <- lapply(seq_len(reps), loop)
        }
    }
    # simplify method
    plot_richness <- do.call(rbind, res)
    pcutils::del_ps("foreach", "doSNOW")

    result <- plot_richness %>%
        dplyr::group_by(., rare, sample) %>%
        dplyr::summarise(., mean(alpha), sd(alpha)) %>%
        dplyr::select(everything(), alpha = 3, sd = 4)
    result$sample <- result$sample %>% factor(., levels = colnames(otutab))

    attributes(result)$type <- "species"
    attributes(result)$mode <- mode
    attributes(result)$method <- method
    class(result) <- c("rare_res", class(result))
    return(result)
}

#' Calculate Abundance-occupancy_relationship
#'
#' @param ... add
#' @param otutab object, such as data.frame or pc_otu
#'
#' @export
#' @return AOR
#' @references
#' Barberán, A., Bates, S. T., Casamayor, E. & Fierer, N. (2012) Using network analysis to explore co-occurrence patterns in soil microbial communities.
#' @examples
#' data(otutab, package = "pcutils")
#' aor(otutab) -> AOR
#' plot(AOR)
aor <- function(otutab, ...) {
    UseMethod("aor")
}

#' @param otutab otutab
#' @param top_r percentage of top relative abundance
#' @param ocup_n percentage of top occupied
#' @param special_n how many occupancy define as specialists
#'
#' @method aor data.frame
#' @rdname aor
#' @exportS3Method
aor.data.frame <- function(otutab, top_r = 0.7,
                           ocup_n = ceiling(0.8 * ncol(otutab)),
                           special_n = ceiling(0.1 * ncol(otutab)), ...) {
    abundance <- relative_abund <- cum <- occupancy <- NULL
    lib_ps("dplyr", library = FALSE)
    otutab %>%
        dplyr::transmute(
            abundance = rowSums(.), relative_abund = abundance / sum(abundance),
            occupancy = rowSums(. > 0)
        ) %>%
        dplyr::arrange(-abundance) -> AOR
    AOR %>% dplyr::mutate(log_abund = log10(relative_abund)) -> AOR

    y <- AOR$occupancy
    x <- AOR$log_abund
    tmp <- try({
        fit <- stats::nls(y ~ stats::SSlogis(x, Asym, xmid, scal), data = data.frame(x, y))
        # summary(fit)
        AOR$fit_y <- stats::predict(fit, newdata = data.frame(x = x))
    })
    if ("try-error" %in% class(tmp)) {
        warning("fit failed")
        AOR$fit_y <- 0
    }

    # find core otus
    AOR %>% dplyr::mutate(cum = cumsum(relative_abund), core = ifelse((cum < top_r) & (occupancy > ocup_n), "core", "others")) -> AOR
    # habitat generalists and specialists
    AOR %>% dplyr::mutate(habitat = ifelse((cum < top_r) & (occupancy > ocup_n), "generalists",
        ifelse((cum < top_r) & (occupancy < special_n), "specialists", "others")
    )) -> AOR

    class(AOR) <- c("AOR", "data.frame")
    return(AOR)
}


#' Plot a AOR
#'
#' @param x AOR object
#' @param ... add
#'
#' @import scales
#' @rdname aor
#' @exportS3Method
#' @method plot AOR
#' @return ggplot
plot.AOR <- function(x, ...) {
    AOR <- x
    log_abund <- occupancy <- fit_y <- core <- habitat <- NULL
    lib_ps("scales", library = FALSE)

    p1 <- ggplot() +
        geom_point(data = AOR, aes(x = 10^log_abund, y = occupancy), col = "#79A5C9", alpha = 0.6, size = 2) +
        geom_line(data = AOR, aes(x = 10^log_abund, y = fit_y), col = "red", size = 1.5) +
        scale_x_log10(labels = scales::label_log()) +
        theme_bw() +
        labs(x = "mean relative abundnace (%)", y = "occupancy")
    # scale_x_continuous()

    p2 <- ggplot() +
        geom_point(data = AOR, aes(x = 10^log_abund, y = occupancy, col = core), alpha = 0.6, size = 2) +
        scale_x_log10(labels = scales::label_log()) +
        scale_color_manual(values = c("others" = "grey", "core" = "#79A5C9")) +
        theme_bw() +
        labs(x = "mean relative abundnace (%)", y = "occupancy")

    p3 <- ggplot() +
        geom_point(data = AOR, aes(x = 10^log_abund, y = occupancy, col = habitat), alpha = 0.6, size = 2) +
        scale_x_log10(labels = scales::label_log()) +
        scale_color_manual(values = c("others" = "grey", "generalists" = "#79A5C9", "specialists" = "#EE79A5")) +
        theme_bw() +
        labs(x = "mean relative abundnace (%)", y = "occupancy")

    return(list(p1, p2, p3))
}
