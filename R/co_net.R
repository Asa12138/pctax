#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param parallel open parallel mode?(default:F)
#' @param threads parallel mode threads
#' @param filename the prefix of saved files
#'
#' @importFrom  ggcor correlate
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' c_net_cal(totu) -> corr
#' sample(t(otutab) %>% as.data.frame(), 50) -> totu2
#' colnames(totu2) <- paste0("tt", 1:50)
#' c_net_cal(totu, totu2) -> corr2
c_net_cal <- function(totu, totu2 = NULL, filename = "occor",parallel=F,threads=4) {
  if (is.null(totu2)) {
    #corr <- suppressWarnings(ggcor::correlate(totu, method = "spearman", cor.test = TRUE, p.adjust = T, p.adjust.method = "BH"))
    if(parallel)corr<-par_cor(totu,threads = threads,method = "spearman",p.adjust.methods=NULL)
    else corr <- suppressWarnings(ggcor::correlate(totu, method = "spearman", cor.test = TRUE, p.adjust = F))
  }
  else {
    #corr <- suppressWarnings(ggcor::correlate(totu, totu2, method = "spearman", cor.test = TRUE, p.adjust = T, p.adjust.method = "BH"))
    corr <- suppressWarnings(ggcor::correlate(totu, totu2, method = "spearman", cor.test = TRUE, p.adjust = F))
  }

  r <- round(corr$r, 5)
  p.value <- round(corr$p.value, 5)
  rownames(p.value) <- rownames(r)
  colnames(p.value) <- colnames(r)
  # save the correlation result
  write.csv(r, paste0(filename, "_r.csv"))
  write.csv(p.value, paste0(filename, "_p.csv"))
  return(list(r = r, p.value = p.value))
}

#' Parallel calculate spearman correlation for one t(otutab)
#'
#' @param totu t(otutab)
#' @param threads parallel mode threads
#' @param method spearman or not
#' @param p.adjust.methods how to adjust p-value (default:NULL, e.g.fdr)
#'
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' par_cor(totu,4) -> corr
#'
par_cor <- function(totu,threads,method = "spearman",p.adjust.methods=NULL){
  lib_ps("foreach")
  lib_ps("doSNOW")
  if(method=="spearman")totu <- apply(totu,2,rank)
  r <- \(rx,ry){
    n <- length(rx)
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    return(c(r,p))
  }
  nc <- ncol(totu)

  pb <- txtProgressBar(max =nc, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
    corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=nc)
    for(j in 1:nc) {
      if(j > i) corr1[,j] <- r(totu[,i],totu[,j])
    }
    corr1
  }
  stopCluster(cl)

  simplify2array(corr)->corr

  rr <- corr[1,,]
  rr <- rr+t(rr)
  diag(rr) <- 1
  pp <- corr[2,,]
  if(!is.null(p.adjust.methods)){
    lp <- lower.tri(pp)
    pa <- pp[lp]
    pa <- p.adjust(pa, p.adjust.methods)
    pp[lower.tri(pp, diag = FALSE)] <- pa
  }
  pp <- pp+t(pp)
  rownames(pp) <- colnames(totu)
  colnames(pp) <- colnames(totu)
  rownames(rr) <- colnames(totu)
  colnames(rr) <- colnames(totu)

  return(list(r = rr,p.value = pp))
}

#' Import corr from .csv file
#'
#' @param filename
#'
#' @return a list with 7 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'BH'}
#' @export
#'
#' @examples
#' \dontrun{
#' # input_corr("occor")->corr
#' }
input_corr <- function(filename) {
  r <- read.csv(paste0(filename, "_r.csv"), row.names = 1,check.names = F)
  p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1,check.names = F)
  return(list(r = r, p.value = p.value))
}

mmscale<-function(x,min_s=0,max_s=1){
  min_s+(x-min(x))/(max(x)-min(x))*(max_s-min_s)
}

#' Construct a network for table corr
#'
#' @param r_thres r_threshold (default:>0.6)
#' @param p_thres p_threshold (default:<0.05)
#' @param del_single should delete single vertexes?
#' @param corr result from c_net_cal()
#'

#' @return an igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' c_net_cal(totu) -> corr
#' c_net_build(corr) -> co_net
#'
c_net_build <- function(corr, r_thres = 0.6, p_thres = 0.05, del_single = T) {
  suppressMessages(lib_ps("igraph"))
  # set thresholds to construct
  occor.r <- corr$r
  occor.p <- corr$p.value
  occor.r[occor.p > p_thres | abs(occor.r) < r_thres] <- 0
  corr <- occor.r
  # make igraph
  flag <- \(corr){
    if (!nrow(corr) == ncol(corr)) {
      return(FALSE)
    }
    if (!all(rownames(corr) == colnames(corr))) {
      return(FALSE)
    }
    return(TRUE)
  }
  if (flag(corr)) {
    go <- graph_from_adjacency_matrix(as.matrix(corr), mode = "undirected", weighted = T, diag = F)
  } else {
    go <- graph_from_incidence_matrix(as.matrix(corr), directed = F, weighted = T)
  }
  # delete single vertexes?
  if (del_single) go <- delete.vertices(go, V(go)[igraph::degree(go) == 0])
  # set vertex attributes
  # set vertices shape
  V(go)$group <- ifelse(V(go)$name %in% rownames(corr), "group1", "group2")
  V(go)$shape <- ifelse(V(go)$group == "group1", "circle", "square")
  V(go)$color <- ifelse(V(go)$group == "group1", "#F095AA", "#96BEEB")
  V(go)$class <- ifelse(V(go)$group == "group1", "group1", "group2")
  V(go)$size <- 5
  # abs edges weight
  go.weight <- E(go)$weight
  E(go)$cor <- go.weight
  E(go)$weight <- abs(go.weight)
  # set edges color
  E.color <- factor(ifelse(go.weight > 0, "#E85D5D", "#48A4F0"), levels = c("#E85D5D", "#48A4F0"))
  E(go)$color <- as.character(E.color)
  # set edges width
  E(go)$width <- 1.5 * E(go)$weight
  return(go)
}

#' Use dataframe to annotate vertexes of a igraph
#'
#' @param go a igraph object
#' @param anno_tab a dataframe using to annotate(with rowname or a name column)

#' @return a annotated igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab")
#' anno_vertex(co_net, taxonomy)
anno_vertex <- function(go, anno_tab) {
  lib_ps("igraph")
  vertex.attributes(go) %>% as.data.frame() -> v_atr
  if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
  v_atr <- left_join(v_atr, anno_tab, by = "name", suffix = c(".x", ""))
  v_atr %>% select(!ends_with(".x")) -> v_atr
  as.list(v_atr) -> vertex.attributes(go)
  return(go)
}


#' Set basic attributes from totu table
#'
#' @param go a igraph object
#' @param totu t(otu table)
#' @param anno_col a dataframe using to annotate(with rowname or a name column)
#' @param totu2 t(otu table)
#' @param anno_col2 a dataframe using to annotate(with rowname or a name column)
#'
#' @import  RColorBrewer
#' @return a igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' c_net_cal(totu) -> corr
#' sample(t(otutab) %>% as.data.frame(), 50) -> totu2
#' colnames(totu2) <- paste0("tt", 1:50)
#'
#' data("c_net")
#' co_net <- c_net_set(co_net, t(otutab), taxonomy %>% select(Phylum))
#'
#' c_net_set(co_net2, totu, taxonomy %>% select(Phylum), totu2) -> co_net2
#'
c_net_set <- function(go, totu = NULL, anno_col = NULL, totu2 = NULL, anno_col2 = NULL) {
  lib_ps("igraph")
  # set vertex attributes
  # set vertices shape
  V(go)$group <- ifelse(V(go)$name %in% colnames(totu), "group1", "group2")
  V(go)$shape <- ifelse(V(go)$group == "group1", "circle", "square")

  # set vertices size
  if (!is.null(totu)) {
    otu_pro <- totu %>%
      colSums() %>%
      data.frame(name = names(.), abundant = ., size = mmscale(., 5, 10))
    if (!is.null(totu2)) {
      otu_pro2 <- totu2 %>%
        colSums() %>%
        data.frame(name = names(.), abundant = ., size = mmscale(., 5, 10))
      otu_pro <- rbind(otu_pro, otu_pro2)
    }
    anno_vertex(go, otu_pro) -> go
  }
  # set vertices color
  if (!is.null(anno_col)) {
    # stopifnot(ncol(anno_col)==1)
    anno_col <- data.frame(name = rownames(anno_col), class = anno_col[, 1])
    if (!is.null(anno_col2)) {
      # stopifnot(ncol(anno_col2)==1)
      anno_col <- data.frame(name = rownames(anno_col2), class = anno_col2[, 1])
      anno_col <- rbind(anno_col, anno_col2)
    }
    if (is.numeric(anno_col$class)) stop("anno is numeric!!!")
    go.col <- droplevels(as.factor(anno_col$class))
    levels(go.col) <- colorRampPalette(brewer.pal(12, "Set3"))(nlevels(go.col))

    anno_col$color <- as.character(go.col)
    anno_vertex(go, anno_col) -> go
  }
  # fill NAs
  V(go)$color <- ifelse(is.na(V(go)$color), "#96BEEB", V(go)$color)
  V(go)$class <- ifelse(is.na(V(go)$class), "group2", V(go)$class)
  V(go)$size <- ifelse(is.na(V(go)$size), 5, V(go)$size)
  # set edge intra-inter
  same <- \(x){
    return(all(x == x[1]))
  }
  E(go)$class <- get.edgelist(go) %>% apply(., 1, \(x)ifelse(same(x %in% colnames(totu)), "intra", "inter"))
  if (length(E(go)$class %>% unique()) > 1) {
    E(go)$lty <- ifelse(E(go)$class == "intra", 1, 5)
  } else {
    E(go)$lty <- 1
  }

  return(go)
}

#' Plot a cor-network
#'
#' @param go a igraph object
#' @param coors the coordinates you saved
#' @param ... additional parameters for plot.igraph()

#' @return a network plot
#' @export
#'
#' @examples
#' data("c_net")
#' c_net_plot(co_net)
#' c_net_plot(co_net2)
#' c_net_plot(co_net3)
c_net_plot <- function(go, coors = NULL, ...) {
  lib_ps("igraph")
  # mode1，普通的单网络
  if (is.null(coors)) coors <- layout.fruchterman.reingold(go)
  set.seed(123)
  vertex.attributes(go) %>% as.data.frame() -> tmp_v
  plot(go, ...,
    main = "Co-occurrence network", layout = coors,
    vertex.label.font = 1, vertex.label.color = "black",
    vertex.label.cex = 0.05 * V(go)$size,
    vertex.label = ifelse(V(go)$size > 8, V(go)$name, NA),
    edge.curved = TRUE, margin = c(0, 0, 0, 0)
  )
  table(E(go)$cor>0)->numofe

  nfe<-\(x) ifelse(is.na(x),0,x)
  legend(1.2, 1, cex = 0.7, legend = c(paste0("+: ",nfe(numofe["TRUE"])), paste0("-: ",nfe(numofe["FALSE"]))),
         col = c("#E85D5D", "#48A4F0"), bty = "n", title = "Correlation", lty = 1)
  f1 <- length(tmp_v$group %>% unique()) > 1
  f2 <- length(tmp_v$class %>% unique()) > 1
  f3 <- length(E(go)$class %>% unique()) > 1

  if (!f1 && f2) {
    legend(-2, 1,
      cex = 0.7, legend = tmp_v$class %>% unique(),
      col = "black", pt.bg = tmp_v$color %>% unique(), bty = "n", title = "class", pch = 21
    )
  }
  # mode2，两个表交互节点
  if (f1 & !f3) {
    legend(-2, 1,
      cex = 0.7, legend = (filter(tmp_v, group == "group1"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", title = "Tab1", pch = 21
    )
    legend(-2, -.5,
      cex = 0.7, legend = (filter(tmp_v, group == "group2"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", title = "Tab2", pch = 22
    )
  }
  # mode3，两个网络本身加上网络间交互
  if (f1 & f3) {
    legend(-2, 1,
      cex = 0.7, legend = (filter(tmp_v, group == "group1"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", title = "Tab1", pch = 21
    )
    legend(-2, -0.5,
      cex = 0.7, legend = (filter(tmp_v, group == "group2"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", title = "Tab2", pch = 22
    )
    table(E(go)$class=="intra")->numofec
    legend(1.2, .7, cex = 0.7, legend = c("intra", "inter"),legend = c(paste0("intra: ",nfe(numofec["TRUE"])), paste0("inter: ",nfe(numofec["FALSE"]))),
           col = "black", bty = "n", title = "", lty = c(1, 5))
  }
}

#' Calculate natural_connectivity
#'
#' @param p a igraph object

#' @return natural_connectivity (numeric)
#' @export
#'
#' @examples
#' make_ring(10) %>% nc()
nc <- function(p) {
  lib_ps("igraph")
  # 获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
  adj_matrix <- as.matrix(as_adj(p, sparse = F))
  adj_matrix[abs(adj_matrix) != 0] <- 1

  # 矩阵的特征分解，获取特征值 λ
  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)

  # 计算“平均特征根”，获得自然连通度
  lambda_sum <- 0
  N <- length(lambda)
  for (i in 1:N) lambda_sum <- lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum / N, base = exp(1))
  lambda_average
}


#' Calculate all indexs of a network
#'
#' @param p a igraph object
#' @param mode calculate what? c("v", "e", "n", "all")

#' @return a 3-elements list
#' \item{n_index}{indexs of the whole network}
#' \item{v_index}{indexs of each vertex}
#' \item{e_index}{indexs of each edge}
#' @export
#'
#' @examples
#' make_graph("Walther") %>% net_par()
net_par <- function(p, mode = c("v", "e", "n", "all")) {
  lib_ps("igraph")
  stopifnot(is_igraph(p))
  if ("all" %in% mode) mode <- c("v", "e", "n")

  n_index <- NULL
  v_index <- NULL
  e_index <- NULL
  up <- p
  if(!is.null(edge.attributes(up)$weight))up<-delete_edge_attr(up,"weight")
  if ("n" %in% mode) {
    # Calculate Network Parameters
    nnodes <- length(V(p)) # number of nodes
    nedges <- length(E(p)) # number od edges
    edensity <- edge_density(p) # density of network, connectance
    aplength <- average.path.length(up) # Average path length
    #w_aplength <- ifelse(is.null(E(p)$weight), aplength, average.path.length(p)) # weighted Average path length
    degree <- mean(igraph::degree(p)) # Average degree
    w_degree <- ifelse(is.null(E(p)$weight), degree, sum(E(p)$weight) / length(V(p))) # weighted degree
    diameter <- diameter(up) # network diameter
    #w_diameter <- ifelse(is.null(E(p)$weight), diameter, diameter(p))
    clusteringC <- transitivity(p) # Clustering coefficient
    cb <- centralization.betweenness(p)$centralization # 介数中心性(Betweenness centralization)
    n_conn <- nc(p) # 自然连通度

    n_index <- data.frame(
      nnodes, nedges, edensity, aplength, degree, w_degree,diameter,clusteringC,cb,n_conn)

if(F){
    # mean_dist=mean_distance(p)#平均距离
    # w_mean_dist=ifelse(is.null(E(p)$weight),mean_dist,mean_distance(p))
    v_conn <- vertex.connectivity(p) # 点连通度
    e_conn <- edge.connectivity(p) # 边连通度
    components <- count_components(p) # 连通分量的数目
    cc <- centralization.closeness(p)$centralization # 紧密中心性
    cd <- centralization.degree(p)$centralization # 度中心性(Degree centralization)
    ce <- centralization.evcent(p)$centralization # 特征向量中心性

    n_index <- data.frame(
      nnodes, nedges, edensity, aplength, degree, w_degree,
      diameter, clusteringC, cb, cc, cd, ce,
      components, n_conn
    )
}

    n_index <- apply(n_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(graph_attr(p) %>% unlist() %>% is.null())) n_index <- data.frame(graph_attr(p), n_index)
  }
  if ("v" %in% mode) {
    # Calculate Vertices Parameters
    v_index <- data.frame(
      degree = degree(p),
      clusteringC = transitivity(p, type = "local"), # 局部聚类系数
      betweenness = betweenness(p), # 节点介数
      eccentricity = eccentricity(p),
      page_rank = page.rank(p)$vector
    )
    v_index <- apply(v_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(vertex_attr(p) %>% unlist() %>% is.null())) v_index <- data.frame(vertex_attr(p), v_index)
  }
  if ("e" %in% mode) {
    # Calculate Edges Parameters
    e_index <- data.frame(
      edge_attr(p)
    )
    # if(!(edge_attr(p)%>%unlist()%>%is.null()))e_index=data.frame(edge_attr(p),e_index)
  }
  return(list(n_index = n_index, v_index = v_index, e_index = e_index))
}


#' Detect the modules
#'
#' @param go a igraph object
#' @param method cluster_method:("cluster_walktrap","cluster_edge_betweenness","cluster_fast_greedy","cluster_spinglass")
#' @return a 3-elements list
#' \item{go}{a igraph object}
#' \item{wc}{a community object}
#' \item{indexs}{net_par result}
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
modu_dect <- function(go, method = "cluster_fast_greedy") {
  lib_ps("igraph")
  stopifnot(is_igraph(go))
  ms <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  if (!method %in% ms) stop("method should be one of ", paste(ms, collapse = ","))

  if (method == "cluster_walktrap") {
    wc <- igraph::cluster_walktrap(go, weights = abs(igraph::E(go)$weight))
  }
  if (method == "cluster_edge_betweenness") {
    wc <- igraph::cluster_edge_betweenness(go, weights = abs(igraph::E(go)$weight))
  }
  if (method == "cluster_fast_greedy") {
    wc <- igraph::cluster_fast_greedy(go, weights = abs(igraph::E(go)$weight), )
  }
  if (method == "cluster_spinglass") {
    wc <- igraph::cluster_spinglass(go, weights = abs(igraph::E(go)$weight))
  }
  V(go)$module <- membership(wc)
  indexs <- net_par(go)
  indexs$n_index$modularity <- modularity(wc)
  # =====模块内模块间分析
  within <- within_module_deg_z_score(go)
  pc <- part_coeff(go)
  indexs$v_index <- data.frame(indexs$v_index, within, pc, check.rows = T)
  {
    taxa.roles <- indexs$v_index
    taxa.roles[
      which(taxa.roles$Zi < 2.5 & taxa.roles$Pi < 0.62),
      "roles"] <- "Peripherals"
    taxa.roles[
      which(taxa.roles$Zi < 2.5 & taxa.roles$Pi >= 0.62),
      "roles"] <- "Connectors"
    taxa.roles[
      which(taxa.roles$Zi >= 2.5 & taxa.roles$Pi < 0.62),
      "roles"] <- "Module hubs"
    taxa.roles[
      which(taxa.roles$Zi >= 2.5 & taxa.roles$Pi >= 0.62),
      "roles"] <- "Network hubs"
    indexs$v_index <- taxa.roles
  }
  return(list(go = go, wc = wc, indexs = indexs))
}

#' Plot a community(modularity)
#'
#' @param wc a community object or a list contain wc and go
#' @param go a igraph object
#' @param n_modu >n_modu modules will be ploted
#' @param ... additional
#' @return a plot
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' modu_plot(co_net_modu, n_modu = 5)
modu_plot <- function(wc, go = NULL, coors=NULL,n_modu = 0, ...) {
  if ("list" %in% class(wc)) go <- wc$go;wc <- wc$wc
  stopifnot("communities" %in% class(wc))
  stopifnot("igraph" %in% class(go))
  as_group <- \(x){
    x <- list(membership = x)
    vids <- names(x$membership)
    modus <- tapply(vids, x$membership, simplify = FALSE, function(x) x)
    return(modus)
  }
  set.seed(123)
  if (is.null(coors)) coors <- layout.fruchterman.reingold(go)
  if (n_modu) {
    # 筛选成员大于4的模块
    membership(wc)[membership(wc) %in% (which((membership(wc) %>% table()) > n_modu))] %>% as_group() -> new_modu
    plot(go,...,
      main = "Modularity network", layout = coors,
      mark.groups = new_modu, vertex.color = V(go)$module,
      vertex.label.font = 2, vertex.label.color = "black",
      vertex.label.cex = .05* V(go)$size, vertex.size = V(go)$size,
      vertex.label = ifelse(V(go)$size > 8, V(go)$name, NA),
      edge.lty = 1, edge.curved = TRUE, margin = c(0, 0, 0, 0)
    )
    legend(1, 1, cex = .7, legend = names(new_modu), fill = rainbow(length(new_modu),
      alpha = 0.3
    ), bty = "n", title = "Modules")
  }
  else {
    plot(wc, go,...,
      main = "Modularity network", layout = coors,
      vertex.label.font = 2, vertex.label.color = "black",
      vertex.label.cex = .05* V(go)$size, vertex.size = V(go)$size,
      vertex.label = ifelse(V(go)$size > 8, V(go)$name, NA),
      edge.lty = 1, edge.curved = TRUE, margin = c(0, 0, 0, 0)
    )
    legend(1.2, 1,
      cex = .7, legend = names(communities(wc)), fill = rainbow(length(communities(wc)), alpha = 0.3),
      bty = "n", title = "Modules"
    )
  }
}


# 计算模块内连接度的z-scores函数
within_module_deg_z_score <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, names = TRUE, attr = "weight")
    } else {
      A <- as_adj(g, sparse = FALSE, names = TRUE)
    }
  }
  memb <- vertex_attr(g, "module")
  N <- max(memb)
  nS <- tabulate(memb)
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  names(z) <- names(Ki) <- rownames(A)
  for (S in seq_len(N)) {
    x <- rowSums(as.matrix(A[memb == S, memb == S]))
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S] - 1))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  z[is.nan(z)] <- 0
  Zi <- z
  df <- data.frame(Ki, Zi, row.names = names(Ki))
  return(df)
}

# 计算节点的参与系数函数
part_coeff <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, attr = "weight")
    } else {
      A <- as_adj(g, sparse = FALSE)
    }
  }
  memb <- vertex_attr(g, "module")
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  Pi <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  names(Pi) <- rownames(A)
  Pi <- data.frame(Pi)
  return(Pi)
}

#' Zi-Pi plot of vertexes
#'
#' @param indexs indexs in modu_dect() result
#'
#' @return a ggplot object
#' @export
#' @import ggpubr ggrepel
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' zp_plot(co_net_modu$indexs)
zp_plot <- function(indexs) {
  taxa.roles <- indexs$v_index
  x1 <- c(0, 0.62, 0, 0.62)
  x2 <- c(0.62, 1, 0.62, 1)
  y1 <- c(-Inf, 2.5, 2.5, -Inf)
  y2 <- c(2.5, Inf, Inf, 2.5)
  lab <- c("peripheral", "Connectors", "Module hubs", "Network hubs")
  CPCOLS <- c("#FC5D6096", "#9BC799B9", "#94CCF2AC", "#FCF6EFFC")
  p <- ggplot() +
    geom_rect(data = NULL, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab)) +
    guides(fill = guide_legend(title = "Topological roles")) +
    scale_fill_manual(values = CPCOLS) +
    geom_point(data = taxa.roles, aes(x = Pi, y = Zi, color = factor(module))) +
    ggrepel::geom_text_repel(
      data = filter(taxa.roles, roles != "Peripherals"),
      aes(x = Pi, y = Zi, label = name)
    ) +
    ggpubr::theme_pubr(legend = "right") +
    guides(colour = "none") +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Participation Coefficient") +
    ylab("Within-module connectivity z-score")
  return(p)
}


#' Robustness_test for a network
#'
#' @param go a igraph object
#' @param partial how much percent vertexes be removed
#'
#' @return dataframe(robustness class)
#' @export
#'
#' @examples
#' data("c_net")
#' robustness_test(co_net,threads=1)->robust_res
#' class(robust_res)
#' plot(robust_res)

robustness_test <- function(go, partial = 0.7,step=1,reps=9,threads=4) {
  cal_del<-\(go,partial,step,rep){
    nodes <- length(V(go))
    floor(nodes * partial) -> del_i
    del_i_indexs <- data.frame()
  for (i in seq(0,del_i,step)) {
    # 在网络中随机移除 i 个节点
    remove_node <- sample(1:nodes, i)
    dp <- delete.vertices(go, remove_node)
    # 计算指标
    del_i_indexs <- rbind(
      del_i_indexs,
      data.frame(net_par(dp,mode = "n")$n_index, i = i)
    )
  }
    return(data.frame(del_i_indexs,"rep"=rep))
  }

  if(threads==1){
    lapply(1:reps, \(i)cal_del(go,partial,step,i))->a
    do.call(rbind,a)->del_i_indexs
  }
  else if(threads>1){
    #parallel
    lib_ps("foreach","doSNOW")
    pb <- txtProgressBar(max =reps, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    cl <- makeCluster(threads)
    registerDoSNOW(cl)
    del_i_indexs<- foreach (i = 1:reps,.combine = "rbind",
                                         .options.snow = opts,
                                         .packages = c("igraph","pcutils")) %dopar%{
                                           cal_del(go,partial,step,i)
                                         }
    stopCluster(cl)
    }
  class(del_i_indexs)<-c("robustness",class(del_i_indexs))
  return(del_i_indexs)
}


#' Plot robustness
#'
#' @param robust_res robustness_test() result (robustness class)
#' @param indexs indexs selected to show
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return a ggplot
#' @method plot robustness
#' @export
#' @rdname robustness_test
plot.robustness<-function(robust_res, indexs = c("n_conn", "aplength", "degree")){
  robust_res %>%
    select(i, all_of(indexs)) %>%
    reshape2::melt(id.var = "i") -> pdat

  pdat%>%group_by(i,variable)%>%summarise(mean=mean(value),sd=sd(value),se=sd/sqrt(length(value)))->sdd
  p <- ggplot(sdd, aes(x=i, y=mean, col = variable)) +
    geom_line()+
    geom_errorbar(data = sdd,aes(ymax=mean+se,ymin=mean-se))+
    #geom_smooth(se = FALSE,method = "loess",formula = 'y ~ x') +
    facet_wrap(. ~ variable, scales = "free")
  return(p)
}

