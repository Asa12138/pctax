# 中文名称==============

#' get all species Latin and Chinese name from the CCTCC database
#'
#' @param download_dir default
#' @param max_requests default 50
#' @param max_id default 30609, try to make sure on the website
#' @param failure_ids failure_ids
#' @param each_verbose each_verbose
#'
#' @return No value
get_all_sp_la_zh_name <- function(download_dir = "~/Documents/", each_verbose = FALSE,
                                  max_requests = 50, max_id = 30609, failure_ids = NULL) {
  lib_ps("httr", library = FALSE)
  lib_ps("jsonlite", library = FALSE)

  # 随机生成 User-Agent 列表
  user_agents <- c(
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
    "Mozilla/5.0 (iPhone; CPU iPhone OS 14_6 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.1.1 Mobile/15E148 Safari/604.1",
    "Mozilla/5.0 (iPad; CPU OS 14_6 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.1.1 Mobile/15E148 Safari/604.1"
  )
  if (!dir.exists(download_dir)) {
    dir.create(download_dir)
  }
  file_dir <- download_dir

  # 初始化存储
  if (file.exists(paste0(file_dir, "all_sp_data.RData"))) load(paste0(file_dir, "all_sp_data.RData"), envir = environment())
  if (!exists("all_sp_data", envir = environment())) all_sp_data <- list()
  if (file.exists(paste0(file_dir, "failure_ids.RData"))) load(paste0(file_dir, "failure_ids.RData"), envir = environment())
  if (!exists("failure_ids", envir = environment())) failure_ids <- c()

  if (is.null(failure_ids)) {
    get_ids <- seq(length(all_sp_data) + 1, max_id)
  } else {
    get_ids <- failure_ids
  }

  # 循环爬取
  for (id in get_ids) {
    if (each_verbose) message(paste0("ID: ", id))
    # 设置 User-Agent 随机值
    user_agent <- sample(user_agents, 1)
    # 构造请求 URL
    api_url <- paste0("http://cctcc.whu.edu.cn/api/dictionary/detail.html?id=", id)

    # 发送请求
    response <- tryCatch(
      {
        httr::GET(api_url, httr::add_headers(`User-Agent` = user_agent))
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(response)) {
      failure_ids <- c(failure_ids, id)
      next
    }

    # 检查响应状态和内容是否为 JSON
    if (httr::status_code(response) == 200) {
      response_text <- httr::content(response, "text", encoding = "UTF-8")
      if (grepl("<HTML><HEAD><TITLE>\u8bbf\u95ee\u7981\u6b62", response_text)) {
        warning("Access is prohibited, wait for 10 minutes...")
        Sys.sleep(600) # 等待 10 分钟
        next
      }

      # 检查 JSON 格式并解析
      if (!grepl("^\\{.*\\}$", response_text)) {
        message(paste0("Non JSON response, skip ID: ", id))
        failure_ids <- c(failure_ids, id)
        next
      }

      # 解析 JSON 数据
      data <- tryCatch(
        {
          jsonlite::fromJSON(response_text)
        },
        error = function(e) {
          message(paste0("Non JSON response, skip ID: ", id))
          failure_ids <- c(failure_ids, id)
          NULL
        }
      )

      if (!is.null(data)) {
        all_sp_data[[id]] <- data$data
      }
    } else {
      failure_ids <- c(failure_ids, id)
    }

    # 随机暂停
    Sys.sleep(runif(1, 0.4, 0.9)) # 随机等待 0.4-0.9 秒

    # 每 50 次保存一次
    if (id %% max_requests == 0) {
      message(paste0("Done ", id, " requests"))
      save(all_sp_data, file = paste0(file_dir, "all_sp_data.RData"))
    }

    # 每 max_requests 次强制暂停
    if (id %% max_requests == 0) {
      Sys.sleep(runif(1, 1, 5)) # 强制等待 1-5 秒
    }
  }

  # 保存失败记录
  if (length(failure_ids) > 0) {
    save(failure_ids, file = paste0(file_dir, "failure_ids.RData"))
  }
  save(all_sp_data, file = paste0(file_dir, "all_sp_data.RData"))
  message("Done!")
}

list_to_dataframe <- function(lst) {
  df <- do.call(rbind, lapply(lst, function(x) {
    x[vapply(x, is.null, logical(length = 1L))] <- NA
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  return(df)
}


save_all_sp_la_zh_name <- function(file = "~/Documents/all_sp_data.RData",
                                   save_file = "~/Documents/R/pctax/pctax/data/all_sp_la_zh_name.rda") {
  all_sp_data <- NULL
  load(file, envir = environment())
  list_to_dataframe(all_sp_data) -> all_sp_df
  all_sp_la_zh_name <- all_sp_df[, c("latin_name", "chinese_name")]
  all_sp_la_zh_name <- dplyr::mutate_all(all_sp_la_zh_name, trimws)
  gsub(
    "\u5206\u7c7b\u540d\u79f0\u53d8\u66f4\u2192",
    "(\u540d\u79f0\u53d8\u66f4)", all_sp_la_zh_name$chinese_name
  ) -> all_sp_la_zh_name$chinese_name
  gsub("\\s+", " ", all_sp_la_zh_name$chinese_name) -> all_sp_la_zh_name$chinese_name
  save(all_sp_la_zh_name, file = save_file)
}

#' Convert taxon names between Chinese and Latin
#'
#' @param input_names input names
#' @param mode conversion mode, "latin_to_chinese" or "chinese_to_latin"
#' @param fuzzy whether to use fuzzy matching, default is FALSE
#'
#' @return character vector of converted names
#' @export
#'
#' @examples
#' convert_taxon_name(c("Escherichia coli", "Clostridioides difficile"))
convert_taxon_name <- function(input_names, mode = "latin_to_chinese", fuzzy = FALSE) {
  # 检查 mode 参数
  if (!mode %in% c("latin_to_chinese", "chinese_to_latin")) {
    stop("mode should be 'latin_to_chinese' or 'chinese_to_latin'")
  }
  all_sp_la_zh_name <- NULL
  data("all_sp_la_zh_name", envir = environment(), package = "pctax")

  mapping_table <- all_sp_la_zh_name
  # 检查映射表是否包含必要列
  if (!("latin_name" %in% colnames(mapping_table)) ||
    !("chinese_name" %in% colnames(mapping_table))) {
    stop("need 'latin_name' and 'chinese_name' columns")
  }

  # 根据模式设置列
  from_col <- if (mode == "latin_to_chinese") "latin_name" else "chinese_name"
  to_col <- if (mode == "latin_to_chinese") "chinese_name" else "latin_name"

  # 转换函数
  convert <- function(input) {
    if (fuzzy) {
      # 模糊匹配
      matches <- mapping_table[grepl(input, mapping_table[[from_col]], ignore.case = TRUE), ]
    } else {
      # 精确匹配
      matches <- mapping_table[mapping_table[[from_col]] == input, ]
    }

    # 返回结果
    if (nrow(matches) == 0) {
      return("") # 未找到匹配
    } else if (nrow(matches) == 1) {
      return(matches[[to_col]])
    } else {
      # 多个匹配，返回第一个
      return(matches[[to_col]][1])
    }
  }

  # 应用转换
  vapply(input_names, convert, FUN.VALUE = character(length = 1L))
}

# phylogenetic tree==============


# 定义目标分类等级及其前缀
rank_prefixes <- c(
  Domain = "d",
  Kingdom = "k",
  Phylum = "p",
  Class = "c",
  Order = "o",
  Family = "f",
  Genus = "g",
  Species = "s"
)

parse_mpa_taxon <- function(tax_string) {
  # 解析函数：把每行的taxonomy分配到正确的列
  parse_taxonomy <- \(tax_string) {
    tax_levels <- strsplit(tax_string, "\\|")[[1]]
    # 初始化空的等级向量
    result <- setNames(rep(NA, length(rank_prefixes)), names(rank_prefixes))

    # 遍历当前层级，把有的放进对应列
    for (level in tax_levels) {
      for (rank in names(rank_prefixes)) {
        prefix <- paste0(rank_prefixes[rank], "__")
        if (startsWith(level, prefix)) {
          result[rank] <- sub(prefix, "", level)
          break
        }
      }
    }
    return(result)
  }

  # 应用到所有行
  tax_parsed <- t(vapply(tax_string, parse_taxonomy, character(length = 8)))
  tax_df <- as.data.frame(tax_parsed, stringsAsFactors = FALSE)
  tax_df
}

#' Load a metaphlan format data.frame
#'
#' @param mpa_df metaphlan format data.frame, rownames are taxon, all value are numeric.
#' @param sum_unidentified logical, whether to sum the unidentified reads to the correspond level.
#'
#' @return a list
#' @export
#'
load_mpa_df <- function(mpa_df, sum_unidentified = TRUE) {
  # 提取分类路径列（第一列）
  parse_mpa_taxon(rownames(mpa_df)) -> tax_df
  # 填充
  pre_tax_table(tax_df, tax_levels = rank_prefixes) -> tax_df2

  tax_ranks <- names(rank_prefixes)

  # 找出每行的最低分类等级
  get_lowest_rank <- \(row) {
    non_na <- which(!is.na(row))
    if (length(non_na) == 0) {
      return(NA)
    }
    return(tax_ranks[max(non_na)])
  }

  make_split_tables <- \(split_tables4){
    for (i in seq_along(tax_ranks)) {
      split_tables4[[i]] <- split_tables4[[i]][, c(tax_ranks[seq_len(i)], sample_cols)]
      rownames(split_tables4[[i]]) <- split_tables4[[i]][, i]
    }
    split_tables4
  }

  # 应用到 tax_df
  lowest_ranks <- apply(tax_df, 1, get_lowest_rank)
  # 添加到总数据中
  tax_df2$LowestRank <- lowest_ranks
  final_df <- cbind(tax_df2, mpa_df)

  LowestRank <- NULL
  # 拆分为7个表格
  split_tables <- lapply(tax_ranks, function(rank) {
    subset(final_df, LowestRank == rank)[, !(names(final_df) %in% "LowestRank")]
  })
  names(split_tables) <- tax_ranks

  sample_cols <- setdiff(colnames(final_df), c(tax_ranks, "LowestRank"))
  split_tables2 <- split_tables

  i <- length(tax_ranks)
  tmp <- hebing(split_tables2[[i]][, sample_cols], group = tax_ranks[i], margin = 1, act = "sum", metadata = split_tables2[[i]])
  tmp2 <- distinct_at(split_tables2[[i]][, seq_along(tax_ranks)], i, .keep_all = TRUE)
  rownames(tmp2) <- tmp2[, i, drop = TRUE]
  split_tables2[[i]] <- cbind(tmp2, tmp[rownames(tmp2), ])

  # 回补上一级分类，un_k__等等
  for (i in seq(length(tax_ranks), 2)) {
    tmp <- hebing(split_tables2[[i]][, sample_cols], group = tax_ranks[i - 1], margin = 1, act = "sum", metadata = split_tables2[[i]])
    tmp2 <- distinct_at(split_tables2[[i]][, seq_along(tax_ranks)], i - 1, .keep_all = TRUE)
    rownames(tmp2) <- tmp2[, i - 1, drop = TRUE]
    tmp3 <- cbind(tmp2, tmp[rownames(tmp2), ])

    huibu_term <- setdiff(rownames(tmp), split_tables2[[i - 1]][, tax_ranks[i - 1]])
    if (length(huibu_term) > 0) {
      split_tables2[[i - 1]] <- rbind(split_tables2[[i - 1]], tmp3[huibu_term, ])
    }
  }

  if (!sum_unidentified) {
    return(make_split_tables(split_tables2))
  }

  # 回补未识别的分类
  split_tables3 <- split_tables2
  for (i in seq(2, length(tax_ranks))) {
    tmp <- hebing(split_tables3[[i]][, sample_cols], group = tax_ranks[i - 1], margin = 1, act = "sum", metadata = split_tables3[[i]])
    rownames(split_tables3[[i - 1]]) <- split_tables3[[i - 1]][, i - 1, drop = TRUE]
    tmp2 <- data.frame(split_tables3[[i - 1]][, sample_cols], check.names = FALSE)

    com_id <- intersect(rownames(tmp2), rownames(tmp))
    diff_df <- tmp2[com_id, ] - tmp[com_id, ]
    tmp3 <- cbind(split_tables3[[i - 1]][com_id, seq_along(tax_ranks)], diff_df)
    tmp3[, i] <- paste0(rank_prefixes[i], "__", tmp3[, i - 1], "_unidentified")

    if (any(tmp3[, sample_cols] > 0)) {
      split_tables3[[i]] <- rbind(split_tables3[[i]], tmp3)
    }
    add_term <- setdiff(rownames(tmp2), rownames(tmp))
    if (length(add_term) > 0) {
      tmp4 <- split_tables3[[i - 1]][add_term, ]
      tmp4[, i] <- paste0(rank_prefixes[i], "__", substr(tmp4[, i - 1], 4, nchar(tmp4[, i - 1])))
      split_tables3[[i]] <- rbind(split_tables3[[i]], tmp4)
    }
  }
  rownames(split_tables3[[i]]) <- split_tables3[[i]][, i, drop = TRUE]
  return(make_split_tables(split_tables3))
}

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
#' @param ignore_pattern An optional regular expression pattern to match tip or node labels for dropping.
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
df2tree <- function(data, edge_df = FALSE, ignore_pattern = NULL) {
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

  if (!is.null(ignore_pattern)) {
    datalist[pcutils::grepl.data.frame(ignore_pattern, datalist)] <- NA
    datalist <- na.omit(datalist)
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


# 注意这个有问题，内部节点的关系会有变化
#' Drop Tips and Update a Phylogenetic Tree
#'
#' This function iteratively removes specified tips (or tips matching a pattern)
#' from a phylogenetic tree without collapsing internal nodes or singleton nodes.
#'
#' @param tr A phylogenetic tree of class \code{phylo}.
#' @param drop_name A character vector of tip or node names to drop.
#'        If missing and \code{pattern} is provided, names matching the pattern will be dropped.
#' @param pattern An optional regular expression pattern to match tip or node labels for dropping.
#'
#' @return A \code{phylo} object with specified tips removed.
#'
#' @examples
#' if (requireNamespace("ape")) {
#'   library(ape)
#'   tr <- rtree(10)
#'   plot(tr)
#'   # Drop tips containing "t1" or "t2" in their label
#'   tr2 <- drop_tips_update(tr, pattern = "t1|t2")
#'   plot(tr2)
#'
#'   # Alternatively, specify tips directly
#'   tr3 <- drop_tips_update(tr, drop_name = c("t3", "t5"))
#'   plot(tr3)
#' }
#' @export
drop_tips_update <- function(tr, drop_name, pattern = NULL) {
  stopifnot(inherits(tr, "phylo"))

  if (missing(drop_name)) {
    if (!is.null(pattern)) {
      drop_name <- c(
        grep(pattern, tr$tip.label, value = TRUE),
        grep(pattern, tr$node.label, value = TRUE)
      )
    } else {
      stop("Either 'drop_name' or 'pattern' must be provided.")
    }
  }

  tr2 <- tr
  while (any(tr2$tip.label %in% drop_name)) {
    tr2 <- ape::drop.tip(tr2, intersect(tr2$tip.label, drop_name),
      trim.internal = FALSE, collapse.singles = FALSE
    )
  }

  return(tr2)
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
#' if (requireNamespace("ape")) {
#'   data(otutab, package = "pcutils")
#'   df2tree1(taxonomy) -> tax_tree
#'   print(tax_tree)
#' }
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
#' @param ignore_pattern An optional regular expression pattern to match tip or node labels for dropping.
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
ann_tree <- function(f_tax, otutab = NULL, level = ncol(f_tax), ignore_pattern = NULL) {
  lib_ps("ggtree", "vctrs", library = FALSE)
  label.x <- label.y <- NULL
  # le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  le <- c("root", colnames(f_tax))
  df2tree(f_tax[, 1:level], ignore_pattern = ignore_pattern) -> tree

  # if (!is.null(drop_name)) {
  #   tree <- drop_tips_update(tree, drop_name = drop_name)
  # }

  tree <- ggtree::fortify(tree)

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
#' if (requireNamespace("ggtree")) {
#'   data(otutab, package = "pcutils")
#'   df2tree(taxonomy[1:50, ]) -> tax_tree
#'   df2tree(taxonomy[51:100, ]) -> tax_tree2
#'   link <- data.frame(from = sample(tax_tree$tip.label, 20), to = sample(tax_tree2$tip.label, 20))
#'   plot_two_tree(tax_tree, tax_tree2, link)
#' }
plot_two_tree <- function(tree1, tree2, edge_df = NULL, tree2_x = 10, filter_link = FALSE,
                          tree1_param = list(), tree2_param = list(),
                          line_param = list(),
                          tree1_tip = FALSE, tip1_param = list(), tree2_tip = FALSE, tip2_param = list(),
                          tree1_highlight = NULL, highlight1_param = list(), highlight1_scale = NULL,
                          tree2_highlight = NULL, highlight2_param = list(), highlight2_scale = ggplot2::scale_fill_hue(na.value = NA)) {
  lib_ps("ggtree", library = FALSE)
  e_type <- group <- isTip <- label <- node <- width <- x <- y <- NULL
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
