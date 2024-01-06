# =====functional gene analysis====
#
#' Gene symbolid transfer to entrezIDs (human gene)
#'
#' @param genes gene symbols e.g:ASGR2
#'
#' @return gene entrezIDs dataframe
#' @export
#'
#' @examples
#' genes <- c(
#'     "ASGR2", "BEST1", "SIGLEC16", "ECRP", "C1QC", "TCN2", "RNASE2",
#'     "DYSF", "C1QB", "FAM20A", "FCGR1A", "CR1", "HP", "VSIG4", "EGR1"
#' )
#' gene2id(genes) -> geneid
gene2id <- function(genes) {
    lib_ps("clusterProfiler", "org.Hs.eg.db", library = F)
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
#'
pre_GEO <- function(my_id, GEO_dir = "GEO_data", file = NULL) {
    # 1.Importing the data
    ## change my_id to be the dataset that you want.
    lib_ps("GEOquery", "Biobase", library = F)

    if (is.null(file)) {
        if (file.exists(paste0(GEO_dir, "/", my_id, "_series_matrix.txt.gz"))) {
            file <- paste0(GEO_dir, "/", my_id, "_series_matrix.txt.gz")
        }
    }

    if (is.null(file)) {
        gse <- GEOquery::getGEO(my_id[1], destdir = GEO_dir)
        gse <- gse[[1]]
    } else {
        gse <- GEOquery::getGEO(filename = file, getGPL = F)
    }

    GPL_version <- gse@annotation

    meta <- Biobase::pData(gse) ## print the sample information
    GPL_data_11 <- Biobase::fData(gse) ## print the gene annotation
    GSE_data_expr <- Biobase::exprs(gse) ## print the expression data

    # 进行注释
    GSE_data_expr <- GSE_data_expr %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(prode_id = "rowname") %>%
        mutate(prode_id = as.character(prode_id))

    if (GPL_version %in% c("GPL570", "GPL571", "GPL96", "GSE22873", "GPL1261", "GPL81")) {
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
    } else if (GPL_version %in% c(GPL1536)) {
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

{
COG_desc <- data.frame(
    Category = c(
        "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING",
        "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING",
        "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING",
        "CELLULAR PROCESSES AND SIGNALING", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM",
        "POORLY CHARACTERIZED", "POORLY CHARACTERIZED"
    ),
    Code = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
    Description = c(
        "Translation, ribosomal structure and biogenesis", "RNA processing and modification", "Transcription", "Replication, recombination and repair", "Chromatin structure and dynamics",
        "Cell cycle control, cell division, chromosome partitioning", "Nuclear structure", "Defense mechanisms", "Signal transduction mechanisms", "Cell wall/membrane/envelope biogenesis",
        "Cell motility", "Cytoskeleton", "Extracellular structures", "Intracellular trafficking, secretion, and vesicular transport", "Posttranslational modification, protein turnover, chaperones",
        "Mobilome: prophages, transposons", "Energy production and conversion", "Carbohydrate transport and metabolism", "Amino acid transport and metabolism", "Nucleotide transport and metabolism",
        "Coenzyme transport and metabolism", "Lipid transport and metabolism", "Inorganic ion transport and metabolism", "Secondary metabolites biosynthesis, transport and catabolism",
        "General function prediction only", "Function unknown"
    ),
    stringsAsFactors = FALSE
)
}

update_ec_info <- function() {
    lib_ps("readxl", library = F)
    ec_node <- readxl::read_excel("~/Desktop/element-cycling.xlsx", sheet = 2) %>%
      arrange_all() %>% as.data.frame()
    ec_link <- readxl::read_excel("~/Desktop/element-cycling.xlsx", sheet = 3) %>%
      arrange_all() %>% as.data.frame()
    ec_gene <- readxl::read_excel("~/Desktop/element-cycling.xlsx", sheet = 1) %>%
      arrange_all() %>% as.data.frame()
    ec_path <- readxl::read_excel("~/Desktop/element-cycling.xlsx", sheet = 4) %>%
      arrange_all() %>% as.data.frame()
    all_ec_info <- list(ec_node = ec_node, ec_link = ec_link, ec_gene = ec_gene, ec_path = ec_path)
    save(all_ec_info, file = "~/Documents/R/pctax/pctax/data/all_ec_info.rda")
}

#' Plot element cycle
#'
#' @param cycle one of c("Carbon cycle","Nitrogen cycle","Phosphorus cycle","Sulfur cycle","Iron cycle")
#' @param anno_df anno_df, columns should contains Gene or KO and Group
#' @param only_anno only show genes in anno_df?
#' @param cell_fill cell fill color
#' @param cell_color cell border color
#' @param chemical_size chemical text size
#' @param chemical_bold chemical text bold
#' @param chemical_color chemical text color
#' @param chemical_label chemical text in geom_label or geom_text?
#' @param reaction_width reaction line width
#' @param reaction_arrow_size reaction arrow size
#' @param reaction_arrow_closed reaction arrow closed?
#' @param gene_or_ko "gene" or "ko"
#' @param gene_size gene text size
#' @param gene_x_offset gene_x_offset
#' @param gene_y_offset gene_y_offset
#' @param gene_label gene text in geom_label or geom_text?
#' @param gene_color gene text color
#' @param gene_bold gene text blod?
#' @param gene_italic gene text italic?
#' @param gene_label_fill gene label fill color
#'
#' @return ggplot
#' @export
#'
#' @examples
#' plot_element_cycle()
plot_element_cycle <- function(cycle = "Nitrogen cycle", anno_df = NULL, only_anno = F, cell_fill = NA, cell_color = "orange",
                               chemical_size = 7, chemical_bold = T, chemical_color = "black", chemical_label = T,
                               reaction_width = 1, reaction_arrow_size = 4, reaction_arrow_closed = T,
                               gene_or_ko = "gene", gene_size = 3, gene_x_offset = 0.3, gene_y_offset = 0.15,
                               gene_label = T, gene_color = NULL, gene_bold = T, gene_italic = T, gene_label_fill = "white") {
    if (file.exists("~/Documents/R/pctax/data/all_ec_info.rda")) {
        load("~/Documents/R/pctax/data/all_ec_info.rda", envir = environment())
    } else {
        data("all_ec_info", package = "pctax", envir = environment())
    }

    ec_node <- all_ec_info$ec_node
    ec_link <- all_ec_info$ec_link
    ec_gene <- all_ec_info$ec_gene
    ec_path <- all_ec_info$ec_path

    ec_node2 <- ec_node %>% filter(Type == cycle)
    ec_node2$Label <- gsub("\\\\", "", ec_node2$Label)

    ec_link2 <- ec_link %>% filter(Type == cycle)
    ec_link2$nudge[is.na(ec_link2$nudge)] <- "0"
    ec_link2$curvature[is.na(ec_link2$curvature)] <- 0

    ec_gene2 <- ec_gene %>% filter(Type == cycle)
    if (gene_or_ko == "ko") ec_gene2$Gene <- ec_gene2$KO
    if (!is.null(anno_df)) {
        if (only_anno) {
            ec_gene2 <- right_join(ec_gene2, anno_df)
        } else {
            ec_gene2 <- left_join(ec_gene2, anno_df)
        }
    }

    # 0.plot cell
    p <- ggplot()

    cells <- tibble::tribble(
        ~Type, ~x1, ~x2, ~y1, ~y2,
        "Carbon cycle", -2.7, 2.7, -2.3, 2.3,
        "Nitrogen cycle", -2.3, 4, -2.3, 2.2,
        "Phosphorus cycle", -2.7, 2.7, -1.5, 2,
        "Sulfur cycle", -2.7, 2.7, -2.3, 2.3,
        "Iron cycle", -2.7, 2.7, -1.2, 1.2
    ) %>%
        as.data.frame() %>%
        column_to_rownames("Type")
    if (cycle %in% rownames(cells)) {
        lib_ps("ggchicklet", library = F)
        p <- p + ggchicklet::geom_rrect(
            data = cells[cycle, ], aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = cell_fill, color = cell_color, size = 2, radius = unit(1, "cm")
        )
    }

    # 1.plot chemicals
    if (chemical_bold) ec_node2$Label <- paste0("bold(", ec_node2$Label, ")")
    if (chemical_label) {
        p1 <- p +
            geom_label(
                data = ec_node2,
                aes(x = X, y = Y, label = paste(Label)), size = chemical_size, color = chemical_color, parse = T
            )
    } else {
        p1 <- p +
            geom_text(
                data = ec_node2,
                aes(x = X, y = Y, label = paste(Label)), size = chemical_size, color = chemical_color, parse = T
            )
    }

    # 2.plot reactions
    p2 <- p1
    for (i in seq_len(nrow(ec_link2))) {
        p2 <- p2 + geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Sub_type),
            linewidth = reaction_width,
            data = ec_link2[i, ], curvature = ec_link2[i, "curvature", drop = T],
            arrow = arrow(length = unit(reaction_arrow_size, "mm"), type = ifelse(reaction_arrow_closed, "closed", "open"))
        )
    }

    # 3.plot genes
    for (i in unique(ec_gene2$Pathway)) {
        if (!i %in% ec_path$Pathway) next
        ec_gene2[which(ec_gene2$Pathway == i), c("X", "Y")] <- pcutils::generate_labels(ec_gene2 %>% filter(Pathway == i) %>% pull(Gene),
            input = unlist(ec_path[which(ec_path$Pathway == i), c("X", "Y")]),
            x_offset = gene_x_offset, y_offset = gene_y_offset
        )
    }
    if (gene_bold & gene_italic) {
        ec_gene2$Gene <- paste0("bolditalic(", ec_gene2$Gene, ")")
    } else if (gene_italic) {
        ec_gene2$Gene <- paste0("italic(", ec_gene2$Gene, ")")
    } else if (gene_bold) ec_gene2$Gene <- paste0("bold(", ec_gene2$Gene, ")")

    if (!is.null(gene_color)) {
        gene_color_param <- list(color = gene_color)
    } else {
        gene_color_param <- list()
    }

    if (!is.null(anno_df)) {
        p3 <- p2 + do.call(geom_label, pcutils::update_param(list(
            data = ec_gene2,
            mapping = aes(x = X, y = Y, label = Gene, color = Sub_type, fill = Group), size = gene_size, parse = T
        ), gene_color_param))
        if (!is.numeric(anno_df$Group)) p3 <- p3 + scale_fill_discrete(na.value = "white")
    } else {
        if (gene_label) {
            p3 <- p2 + do.call(geom_label, pcutils::update_param(list(
                data = ec_gene2,
                mapping = aes(x = X, y = Y, label = Gene, color = Sub_type),
                fill = gene_label_fill, size = gene_size, parse = T, show.legend = F
            ), gene_color_param))
        } else {
            p3 <- p2 + do.call(geom_text, pcutils::update_param(list(
                data = ec_gene2,
                mapping = aes(x = X, y = Y, label = Gene, color = Sub_type),
                size = gene_size, parse = T, show.legend = F
            ), gene_color_param))
        }
    }

    # 4.final plot
    lims <- tibble::tribble(
        ~Type, ~x1, ~x2, ~y1, ~y2,
        "Carbon cycle", -2.8, 2.8, -2.5, 2.5,
        "Nitrogen cycle", -2.5, 4, -3, 2.6,
        "Phosphorus cycle", -2.8, 2.8, -2.5, 2.5,
        "Sulfur cycle", -2.8, 2.8, -2.5, 2.5,
        "Iron cycle", -2.8, 2.8, -2.5, 2
    ) %>%
        as.data.frame() %>%
        column_to_rownames("Type")

    p4 <- p3 + coord_fixed() +
        xlim(lims[cycle, 1], lims[cycle, 2]) + ylim(lims[cycle, 3], lims[cycle, 4]) +
        scale_color_manual(values = pcutils::get_cols(length(unique(ec_link2$Sub_type)), "col3")) +
        pctax_theme + labs(x = NULL, y = NULL) +
        theme(
            legend.position = "right", panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text = element_blank()
        )
    message("recommend ggsave(width = 12,height = 10)")
    p4
}
