#' Batch mediate
#'
#' @param data data.frame with two columns: X (independent variable) and Y (dependent variable).
#' @param mediator_df data.frame with mediators, each column representing a different mediator variable.
#' @param nsims Number of bootstrap simulations for estimating confidence intervals (default is 500).
#' @param conf.level Confidence level for the confidence intervals (default is 0.95).
#'
#' @return data.frame
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' X <- rnorm(n)
#' M1 <- 0.5 * X + rnorm(n)
#' M2 <- 0.3 * X + rnorm(n)
#' M3 <- 0.1 * X + rnorm(n)
#' Y <- 0.3 * X + 0.4 * M1 + 0.2 * M2 + rnorm(n)
#' data <- data.frame(X, Y)
#' mediators <- data.frame(M1, M2, M3)
#' if (requireNamespace("mediation")) {
#'   results <- batch_mediate(data, mediators, nsims = 99)
#'   print(results)
#' }
#' @export
batch_mediate <- function(data, mediator_df, nsims = 500, conf.level = 0.95) {
  # 检查必要的包是否已安装
  lib_ps("mediation", library = FALSE)

  message("Use the first two columns of data as X and Y")
  colnames(data) <- c("X", "Y")

  # 初始化结果数据框
  results <- data.frame(
    Mediator = character(),
    ACME = numeric(), # 平均因果中介效应
    ADE = numeric(), # 平均直接效应
    Total_Effect = numeric(), # 总效应
    Prop_Mediated = numeric(), # 中介比例
    ACME_p = numeric(), # ACME的p值
    ADE_p = numeric(), # ADE的p值
    Total_p = numeric(), # 总效应的p值
    stringsAsFactors = FALSE
  )

  # 循环处理每个中介变量
  for (mediator_name in colnames(mediator_df)) {
    # 合并数据
    full_data <- cbind(data, mediator_df[mediator_name])
    colnames(full_data)[ncol(full_data)] <- "M_temp"

    # 拟合模型
    model.m <- stats::lm(M_temp ~ X, data = full_data[, c("X", "M_temp"), drop = FALSE])
    model.y <- stats::lm(Y ~ ., data = full_data[, c("Y", "X", "M_temp")])

    # 进行中介分析
    med <- mediation::mediate(
      model.m,
      model.y,
      treat = "X",
      mediator = "M_temp",
      boot = TRUE,
      sims = nsims,
      conf.level = conf.level
    )

    # 提取结果
    res_row <- data.frame(
      Mediator = mediator_name,
      ACME = med$d0,
      ADE = med$z0,
      Total_Effect = med$tau.coef,
      Prop_Mediated = med$n0,
      ACME_p = med$d0.p,
      ADE_p = med$z0.p,
      Total_p = med$tau.p,
      stringsAsFactors = FALSE
    )

    # 添加到结果数据框
    results <- rbind(results, res_row)
  }

  # 计算置信区间
  ci_lower <- (1 - conf.level) / 2
  ci_upper <- conf.level + ci_lower

  return(results)
}
