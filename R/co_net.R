#=======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param parallel open parallel mode?(default:F)
#' @param threads parallel mode threads
#' @param filename the prefix of saved files
#'
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' metadata[,12:19] -> env
#' c_net_cal(totu, env) -> corr2
c_net_cal <- function(totu, totu2 = NULL, filename = "occor",threads=4) {
  corr<-par_cor(totu,totu2,threads = threads,method = "spearman",p.adjust.methods=NULL)
  r <- round(corr$r, 5)
  p.value <- round(corr$p.value, 5)
  diag(r)=0
  # save the correlation result
  if(is.character(filename)){
    write.csv(r, paste0(filename, "_r.csv"))
    write.csv(p.value, paste0(filename, "_p.csv"))
  }
  return(list(r = r, p.value = p.value))
}

#' Parallel calculate spearman correlation for one t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param threads parallel mode threads, default 2
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
par_cor <- function(totu, totu2 = NULL,threads=2,method = "spearman",p.adjust.methods=NULL){
  lib_ps("foreach")
  lib_ps("doSNOW")
  if(method=="spearman")totu <- apply(totu,2,rank)
  r_p <- \(rx,ry){
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
  if(is.null(totu2)){
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=nc)
      for(j in 1:nc) {
        if(j > i) corr1[,j] <- r_p(totu[,i],totu[,j])
      }
      corr1
    }
    simplify2array(corr)->corr
    rr <- corr[1,,]
    rr <- rr+t(rr)
    diag(rr) <- 1

    pp <- corr[2,,]
    pp <- pp+t(pp)
    rownames(rr) <-rownames(pp) <- colnames(totu)
    colnames(rr) <-colnames(pp) <- colnames(totu)
  }
  else{
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=ncol(totu2))
      for(j in 1:ncol(totu2)) {
        corr1[,j] <- r_p(totu[,i],totu2[,j])
      }
      corr1
    }
    simplify2array(corr)->corr
    rr <- t(corr[1,,])
    pp <- t(corr[2,,])

    rownames(rr) <-rownames(pp) <- colnames(totu)
    colnames(rr) <-colnames(pp) <- colnames(totu2)
  }
  stopCluster(cl)

  if(!is.null(p.adjust.methods))pp<-p.adjust.table(pp,p.adjust.methods)
  return(list(r = rr,p.value = pp))
}


#' SparCC correlation for a otutab
#'
#' @param totu an t(otutab)
#' @param threads default: 4
#' @param filename saved file name
#'
#' @return list contain a correlation matrix and a bootstrap p_value matrix
#' @export
#'
#' @examples
#' par_sparcc(totu)->sparcc_corr
par_sparcc<-function(totu,filename="sparcc",threads=4){
  #sparcc
  lib_ps("SpiecEasi")

  #执行 sparcc 分析
  set.seed(123)
  totu.sparcc <- sparcc(totu,iter = 10,inner_iter = 5)
  sparcc0 <- totu.sparcc$Cor  #稀疏相关性矩阵
  rownames(sparcc0)=colnames(sparcc0)=colnames(totu)

  #通过 100 次自举抽样获取随机矩阵
  set.seed(123)
  reps = 100
  tmpd=tempdir()

  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  foreach (rep = 1:reps,
           .options.snow = opts,
           .packages = c("SpiecEasi")) %dopar%{
             totu.boot <- sample(totu, replace = TRUE)  #bootstrap
             totu.sparcc_boot <- sparcc(totu.boot,iter = 10,inner_iter = 5)  #sparcc 参数设置和上文保持一致
             sparcc_boot <- totu.sparcc_boot$Cor
             colnames(sparcc_boot)=rownames(sparcc_boot) =colnames(totu.boot)
             write.table(sparcc_boot, paste(tmpd,'/sparcc_boot', rep, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #输出随机值的稀疏相关性矩阵
           }
  stopCluster(cl)
  #基于上述观测值的稀疏相关性矩阵以及 100 次 bootstrap 的结果，获取稀疏相关性的伪 p 值
  p <- sparcc0
  p[p!=0] <- 0
  for (i in 1:reps) {
    p_boot <- read.delim(paste(tmpd,'/sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
    p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
  }
  p <- p / reps

  # save the correlation result
  if(is.character(filename)){
    write.csv(sparcc0, paste0(filename, "_r.csv"))
    write.csv(p, paste0(filename, "_p.csv"))
  }
  #unlink(tmpd,recursive = T,force = T)
  return(list(r = sparcc0, p.value = p))
}


#' p.adjust apply on a table (matrix or data.frame)
#'
#' @param pp table
#' @param method see \code{\link[stats]{p.adjust}}
#'
#' @return matrix
#' @export
#'
#' @examples
#' matrix(abs(rnorm(100,0.01)),10,10)->pp
#' p.adjust.table(pp)
p.adjust.table<-\(pp,method="BH"){
  pp<-as.matrix(pp)
  #self-cor or two-tables cor?
  flag <- \(corr){
    if (!nrow(corr) == ncol(corr)) {
      return(FALSE)
    }
    if (!all(rownames(corr) == colnames(corr))) {
      return(FALSE)
    }
    return(TRUE)
  }
  if(flag(pp)){
    lp <- lower.tri(pp)
    pa <- pp[lp]
    pa <- p.adjust(pa,method)
    pp[lower.tri(pp, diag = FALSE)] <- pa
    pp[upper.tri(pp,diag = F)]<-0
    pp<-pp+t(pp)
  }
  else{
    pp1<-p.adjust(pp,method)
    pp<-matrix(pp1,nrow(pp),ncol(pp),dimnames = list(rownames(pp),colnames(pp)))
  }
  return(pp)
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
  if((max(x)-min(x))==0)return(rep((min_s+max_s)/2,length(x)))
  min_s+(x-min(x))/(max(x)-min(x))*(max_s-min_s)
}

#https://mp.weixin.qq.com/s?__biz=MzIxNzc1Mzk3NQ==&mid=2247486852&idx=1&sn=cf80513c46001f0476f4c1eb4c3d18ee&chksm=97f5bd9ca082348af2652aaa2d43d6039ff1ca0ad0abb72a8e040c6a9002f29aab24e8429a84&cur_album_id=1366880700036726784&scene=190#rd
#偏相关网络分析
if(F){
  library(ppcor)

  y.data <- data.frame(
    hl = c(7, 15, 19, 15, 21, 22, 57, 15, 20, 18),
    disp = c(0.000, 0.964, 0.000, 0.000, 0.921, 0.000, 0.000, 1.006, 0.000, 1.011),
    deg = c(9, 2, 3, 4, 1, 3, 1, 3, 6, 1),
    BC = c(1.78e-02, 1.05e-06, 1.37e-05, 7.18e-03, 0.00e+00, 0.00e+00, 0.00e+00, 4.48e-03, 2.10e-06, 0.00e+00)
  )
  #例如在“y.data”中控制“deg”时，计算 hl 和 disp 的 spearman 偏（半偏）相关系数
  pcor.test(x = y.data$hl, y = y.data$disp, z = y.data$deg, method = 'spearman')  #partial correlation
  spcor.test(x = y.data$hl, y = y.data$disp, z = y.data$deg, method = 'spearman')  #semi-partial correlation

  #例如在“y.data”中控制“deg”和“BC”时，计算 hl 和 disp 的 spearman 偏（半偏）相关系数
  pcor.test(x = y.data$hl, y = y.data$disp, z = y.data[, c('deg', 'BC')], method = 'spearman')  #partial correlation
  spcor.test(x = y.data$hl, y = y.data$disp, z = y.data[, c('deg', 'BC')], method = 'spearman')  #semi-partial correlation
  y.data.pcor <- pcor(t(otutab), method = 'spearman')  #partial correlation，以 spearman 偏相关为例
  y.data.pcor
  y.data.pcor$estimate  #偏相关系数矩阵
  y.data.pcor$p.value  #p 值矩阵



}

#=========2.build======
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
#' t(otutab) -> totu
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
  if (del_single) go <- igraph::delete.vertices(go, V(go)[igraph::degree(go) == 0])
  # set vertex attributes
  # set vertices shape
  V(go)$group <- ifelse(V(go)$name %in% rownames(corr), "group1", "group2")
  V(go)$shape <- ifelse(V(go)$group == "group1", "circle", "square")
  V(go)$color <- ifelse(V(go)$group == "group1", "#A6CEE3","#B15928")
  V(go)$class <- ifelse(V(go)$group == "group1", "group1", "group2")
  V(go)$size <- 5
  # abs edges weight
  go.weight <- E(go)$weight
  E(go)$cor <- go.weight
  E(go)$weight <- abs(go.weight)
  # set edges color
  E(go)$color <- ifelse(go.weight > 0, "#48A4F0", "#E85D5D")
  E(go)$inter<-ifelse(go.weight > 0, "positive", "negative")
  # set edges width
  E(go)$width <- 1.5 * E(go)$weight
  return(go)
}


#' Get RMT threshold for a correlation matrix
#'
#' @param occor.r a correlation matrix
#'
#' @return a r-threshold
#' @export
#'
#' @examples
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' RMT_threshold(corr$r)
RMT_threshold<-function(occor.r){
  #RMT的实现
  lib_ps("RMThreshold")
  #Get threshold by RMT
  nwd=getwd()
  if(!dir.exists("./RMT_temp"))dir.create("./RMT_temp")
  setwd("./RMT_temp")

  RMThreshold::rm.matrix.validation(occor.r,unfold.method = "spline")
  res <- RMThreshold::rm.get.threshold(occor.r)

  setwd(nwd)
  #thre <- (res$sse.chosen + res$p.ks.chosen)/2
  thre <-res$chosen.thresholds
  thre
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
#' @param totu2 t(otu table) if it exist
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
#' metadata[,12:19] -> env
#'
#' data("c_net")
#' co_net <- c_net_set(co_net, t(otutab), taxonomy %>% select(Phylum))
#'
#' c_net_set(co_net2, totu, taxonomy %>% select(Phylum), env) -> co_net2
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
    levels(go.col) <- colorRampPalette(brewer.pal(8, "Set2"))(nlevels(go.col))

    anno_col$color <- as.character(go.col)
    anno_vertex(go, anno_col) -> go
  }
  # fill NAs
  V(go)$color <- ifelse(is.na(V(go)$color), "#B15928", V(go)$color)
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

#========3.layout========
#' Layout coordinates
#'
#' @param go igraph
#' @param niter iterations number
#'
#' @export
#'
c_net_lay<-function(go,niter = 500){
  coors<-layout_with_fr(go,niter = niter,grid = "nogrid")
  coors<-data.frame(name = V(go)$name,X=coors[,1],Y=coors[,2])
  return(coors)
}
if(F){
  #ggClustNet里面的几个布局算法的使用
  #https://blog.csdn.net/qazplm12_3/article/details/125967064
  vertex.attributes(co_net_modu)%>%as.data.frame()%>%select(name,module)%>%rename(ID=name,group=module)->netClu
  netClu$group<-as.factor(netClu$group)
  corr$r[netClu$ID,netClu$ID]->cor
  cor[cor<0.7]=0
  result2 = PolygonMaptreeG(cor = cor,nodeGroup =  netClu)
  result2%>%select(3,1,2)%>%rename(name=elements)->coors
  modu_plot(co_net_modu,coors = coor)
  c_net_plot(co_net_modu,coors = coors)
}


#' Layout with group
#'
#' @param go igraph object
#' @param group group name (default:module)
#' @param zoom1 big network layout size
#' @param zoom2 average sub_network layout size
#' @param layout1 layout1 method, one of
#' {1.numeric: use adjusted circle}
#' {2.a dataframe: rowname is group, two columns are X and Y}
#' {3.function: \code{\link[igraph]{layout_}} default: in_circle()}
#' @param layout2 one of \code{\link[igraph]{layout_}}
#'
#' @return coors
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' g_lay(co_net_modu,group ="module" ,zoom2 = 5,layout2 =nicely())->oridata
#' modu_plot(co_net_modu,coors = oridata)
#' g_lay_nice(co_net_modu,group ="module")->oridata
#' modu_plot(co_net_modu,coors = oridata)
#'
g_lay<-function(go,group="module",zoom1=20,zoom2=3,layout1=in_circle(),layout2=in_circle()){
  stopifnot(is_igraph(go))
  if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
  igraph::as_data_frame(go,what = "vertices")%>%select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  xs = as.data.frame(table(nodeGroup$group))
  r = xs$Freq/10
  scale_f=mmscale(r,0.5,2)
  names(scale_f)=levels(nodeGroup$group)

  big_lay<-\(zoom1,layout1){
    if(is.data.frame(layout1))da=layout1
    if(is.numeric(layout1)){
    #弧度修正
    arg=cumsum(scale_f/sum(scale_f)*360)
    x = rep(0, length(r))
    y = rep(0, length(r))
    for (i in 1:length(r)) {
      x[i] =  sin(arg[i] * 3.14/180)
      y[i] =  cos(arg[i] * 3.14/180)
    }
    da = data.frame(x = x, y = y,row.names = levels(nodeGroup$group))
    }
    else {da = data.frame(layout_(make_ring(length(r)),layout1),
                         row.names = levels(nodeGroup$group))
    }
    #group的中心点
    da[,1]=mmscale(da[,1],-zoom1,zoom1)
    da[,2]=mmscale(da[,2],-zoom1,zoom1)
    return(da)
  }

  da = big_lay(zoom1,layout1)

  #每个group内部的分布
  oridata=data.frame()
  for (i in levels(nodeGroup$group)) {
    nodeGroup%>%filter(group==i)%>%pull(ID)->tmpid
    induced_subgraph(go,tmpid,impl = "copy_and_delete")->tmp_net
    igraph::layout_(tmp_net,layout2)->data

    zoom2*scale_f[i]->r2
    data[,1]=mmscale(data[,1],-r2,r2)
    data[,2]=mmscale(data[,2],-r2,r2)
    data <- data.frame(name =tmpid ,X = data[,1] + da[i, 1], Y = data[,2] +da[i, 2])

    oridata = rbind(oridata, data)
  }
  print("Big layout:")
  print(da)
  return(oridata)
}

#' Layout with group nicely
#'
#' @export
#'
#' @rdname g_lay
g_lay_nice<-function(go,group="module"){
  lib_ps("ggraph")
  stopifnot(is_igraph(go))
  if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
  igraph::as_data_frame(go,what = "vertices")%>%select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  edge = data.frame(from = paste("model_", nodeGroup$group,sep = ""), to = nodeGroup$ID)

  vertices_t <- data.frame(name = unique(c(as.character(edge$from),
                                           as.character(edge$to))))
  vertices_t$size = sample(1:10, nrow(vertices_t), replace = TRUE)

  mygraph <- igraph::graph_from_data_frame(edge, vertices = vertices_t)
  data = ggraph::create_layout(mygraph, layout = "circlepack",weight = size)
  coor = data %>% dplyr::filter(leaf == TRUE) %>% dplyr::select(name,x,y)
  colnames(coor) = c("name", "X", "Y")
  return(coor)
}

#========4.plot========
#' Plot a co-network
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
  else {
    coors<-coors[match(V(go)$name,coors$name),2:3]%>%as.matrix()
  }

  set.seed(123)
  vertex.attributes(go) %>% as.data.frame()-> tmp_v
  #show labels
  tmp_v%>%top_n(10,size)%>%pull(name)->toplabel
  tmp_v$label=ifelse(tmp_v$name%in%toplabel,tmp_v$name,NA)

  plot(go, ...,
    main = "Co-occurrence network", layout = coors,
    vertex.label.font = 1, vertex.label.color = "black",
    vertex.label.cex = 0.05 * V(go)$size,
    vertex.label=tmp_v$label,
    edge.curved = TRUE, margin = c(0, 0, 0, 0)
  )

  table(E(go)$cor>0)->numofe

  nfe<-\(x) ifelse(is.na(x),0,x)
  legend(1.2, 1, cex = 0.7, legend = c(paste0("+: ",nfe(numofe["TRUE"])), paste0("- : ",nfe(numofe["FALSE"]))),
         col = c("#48A4F0","#E85D5D"), bty = "n", title = "Correlation", lty = 1)

  f1 <- length(tmp_v$group %>% unique()) > 1
  f2 <- length(tmp_v$class %>% unique()) > 1
  f3 <- length(E(go)$class %>% unique()) > 1

  if (!f1 && f2) {
    legend(-2, 1,
      cex = 0.7, legend = tmp_v$class %>% unique(),
      col = "black", pt.bg = tmp_v$color %>% unique(), bty = "n", pch = 21
    )
  }
  # mode2，两个表交互节点
  if (f1 & !f3) {
    legend(-2, 1,
      cex = 0.7, legend = (filter(tmp_v, group == "group1"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", pch = 21
    )
    legend(-2, -.5,
      cex = 0.7, legend = (filter(tmp_v, group == "group2"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", pch = 22
    )
  }
  # mode3，两个网络本身加上网络间交互
  if (f1 & f3) {
    legend(-2, 1,
      cex = 0.7, legend = (filter(tmp_v, group == "group1"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", pch = 21
    )
    legend(-2, -0.5,
      cex = 0.7, legend = (filter(tmp_v, group == "group2"))$class %>% unique(),
      col = "black", pt.bg = (filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", pch = 22
    )
    table(E(go)$class=="intra")->numofec
    legend(1.2, .7, cex = 0.7,legend = c(paste0("intra: ",nfe(numofec["TRUE"])), paste0("inter: ",nfe(numofec["FALSE"]))),col = "black", bty = "n", title = "", lty = c(1, 5))
  }
}


#' Input a graphml file exported by Gephi
#'
#' @param file graphml file exported by Gephi
#'
#' @return list contains the igraph object and coordinates
#' @export
#'
input_gephi<-function(file){
  igraph::read.graph(file,format = "graphml")->gephi
  vertex_attr(gephi)%>%as.data.frame()->tmp_v
  #extract coors
  coors<-tmp_v[,c("x","y")]
  coors<-data.frame(name = tmp_v$name,X=coors[,1],Y=coors[,2])
  coors%>%mutate(X=mmscale(X,-40,40),Y=mmscale(Y,-40,40))->coors
  #transfrom color
  rgb2code(tmp_v[,c("r","g","b")])%>%pull(code)->tmp_v$color
  E(gephi)$color=ifelse(E(gephi)$cor > 0, "#48A4F0", "#E85D5D")
  #scale size
  tmp_v$size=mmscale(tmp_v$size,1,5)
  E(gephi)$width=mmscale(E(gephi)$width,0.05,0.2)
  #delete
  tmp_v%>%select(-c("label","x","y","r","g","b","id"))%>%as.list()->vertex.attributes(gephi)
  edge.attributes(gephi)["Edge Label"]=edge.attributes(gephi)["id"]=NULL
  return(list(go=gephi,coors=coors))
}

#' Transfer a igraph object to a ggig
#'
#' @param go igraph
#' @param net_pars net_par() result
#' @param coors coordinates for nodes,columns: name, X, Y
#'
#' @return ggig object
#' @export
#'
#' @examples
#' coors=layout_with_fr(co_net, niter=999,grid="nogrid")
#' coors<-data.frame(name = V(co_net)$name,X=coors[,1],Y=coors[,2])
#' to.ggig(co_net,coors)->ggig
#' plot(ggig)
to.ggig<-function(go,net_pars=NULL,coors=NULL){
  if(is.null(net_pars))net_par(go)->net_pars
  if(is.null(coors)){
    layout.fruchterman.reingold(go)->coors
    coors<-data.frame(name = V(go)$name,X=coors[,1],Y=coors[,2])
  }

  net_pars$v_index%<>%left_join(.,coors,by="name",suffix = c("", ".1"))
  net_pars$e_index%<>%left_join(.,coors,by=c("from"="name"))%>%rename(X1=X,Y1=Y)%>%
    left_join(.,coors,by=c("to"="name"))%>%rename(X2=X,Y2=Y)

  #show labels
  net_pars$v_index%>%top_n(10,size)%>%pull(name)->toplabel
  net_pars$v_index$label=ifelse(net_pars$v_index$name%in%toplabel,net_pars$v_index$name,NA)
  class(net_pars)<-c("ggig","list")
  return(net_pars)
}

#' Plot a ggig
#'
#' @param ggig ggig object
#'
#' @return ggplot
#' @exportS3Method
#'
plot.ggig<-function(ggig){
  lib_ps("ggplot2","ggnewscale")
  p <- ggplot() +
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = inter, linewidth = width),data = ggig$e_index,alpha=0.5) +
    scale_linewidth(range = c(0.2,1.5))+scale_color_manual(values =c(positive="#48A4F0",negative="#E85D5D") )+
    ggnewscale::new_scale_color()+
    geom_point(aes(X, Y,fill = class,size = size),pch = 21, data = ggig$v_index) +
    ggnewscale::new_scale("size")+
    geom_text(aes(X, Y,col = class,size = size,label=label), data = ggig$v_index,show.legend = F)+
    scale_size(range = c(1,4))+
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    coord_fixed(ratio = 1)+
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  p
}


#========5.topological=======
#' Extract sub-nwtwork from the whole network
#'
#' @param a_net the whole network
#' @param otutab otutab, these columns will be extract
#' @param threads threads, default:4
#' @param save_net should save these sub_nets? F or a filename
#'
#' @return a dataframe contains all sub_net parameters
#' @export
#'
#' @examples
#' data(otutab)
#' c_net_cal(t(otutab))->corr
#' corr$p.value<-p.adjust.table(corr$p.value)
#' c_net_build(corr)->a_net
#' extract_sub_net(a_net,otutab,save_net = "testnet")
extract_sub_net<-function(a_net,otutab,threads=4,save_net=F){
  lib_ps("igraph")
  V(a_net)$name->v_name
  reps=ncol(otutab)

  lib_ps("parallel","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  sub_nets <- foreach(i = 1:reps,.packages = c("igraph")) %dopar% {
    rownames(otutab)[otutab[,i]>0]->exist_sp
    induced.subgraph(a_net,which(v_name%in%exist_sp),impl = "copy_and_delete")->spe_sub
  }
  names(sub_nets)<-colnames(otutab)

  sub_net_pars <- foreach(i = 1:reps,.options.snow = opts,
                          .packages = c("igraph","pctax","pcutils"),.combine = "rbind") %dopar% {
                            sub_nets[[i]]->spe_sub
                            net_par(spe_sub,mode = "n")[["n_index"]]->indexs
                            wc <- igraph::cluster_fast_greedy(spe_sub, weights = abs(igraph::E(spe_sub)$weight), )
                            indexs$modularity <- modularity(wc)
                            indexs
                          }

  stopCluster(cl)
  rownames(sub_net_pars)<-colnames(otutab)
  if(is.character(save_net))save(sub_nets,file = paste0(save_net,".rda"))
  sub_net_pars
}

#' Calculate natural_connectivity
#'
#' @param p a igraph object
#' @return natural_connectivity (numeric)
#' @export
#' @references \code{\link[ggClusterNet]{nc}}
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
#' @param go a igraph object
#' @param fast less indexes for faster calculate ?
#' @param mode calculate what? c("v", "e", "n", "all")
#'
#' @return a 3-elements list
#' \item{n_index}{indexs of the whole network}
#' \item{v_index}{indexs of each vertex}
#' \item{e_index}{indexs of each edge}
#' @export
#'
#' @examples
#' make_graph("Walther") %>% net_par()
net_par <- function(go, mode = c("v", "e", "n", "all"),fast=T) {
  pcutils::lib_ps("igraph")
  stopifnot(is_igraph(go))
  if ("all" %in% mode) mode <- c("v", "e", "n")

  n_index <- NULL
  v_index <- NULL
  e_index <- NULL
  #non-weighted network
  up <- go
  if(!is.null(igraph::edge_attr(up)[["weight"]]))up<-delete_edge_attr(up,"weight")
  if ("n" %in% mode) {
    # Calculate Network Parameters
    n_index <- data.frame(
      num_nodes = length(V(go)), # number of nodes
      num_edges = length(E(go)), # number od edges
      edge_density = edge_density(go), # density of network, connectance
      neg_percent=sum(E(go)$cor<0)/length(E(go)), # 负相关比例
      ave_path_len = average.path.length(up), # Average path length
      #w_ave_path_len = ifelse(is.null(E(go)$weight), ave_path_len, average.path.length(go)) # weighted Average path length
      global_efficiency=igraph::global_efficiency(up),
      ave_degree = mean(igraph::degree(go)), # Average degree
      w_ave_degree = ifelse(is.null(E(go)$weight), mean(igraph::degree(go)), sum(E(go)$weight) / length(V(go))), # weighted degree
      diameter = diameter(up), # network diameter
      clusteringC = transitivity(go), # Clustering coefficient
      cen_betweenness = centralization.betweenness(go)$centralization, # 介数中心性(Betweenness centralization)
      nat_connectivity = nc(go) # 自然连通度
    )

    if(!fast){
      # mean_dist=mean_distance(go)#平均距离
      # w_mean_dist=ifelse(is.null(E(go)$weight),mean_dist,mean_distance(go))
      # v_conn= vertex.connectivity(go) # 点连通度
      # e_conn= edge.connectivity(go) # 边连通度
      # components= count_components(go) # 连通分量的数目
      modularity=modularity(cluster_fast_greedy(go))#模块指数
      rand.g <- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
      rand_m=modularity(cluster_fast_greedy(rand.g))
      r_modularity=(modularity-rand_m)/rand_m#相对模块指数

      n_index <- data.frame(
        n_index,
        modularity=modularity,
        r_modularity=r_modularity,
        cen_closeness = centralization.closeness(go)$centralization, # 紧密中心性
        cen_degree = centralization.degree(go)$centralization, # 度中心性(Degree centralization)
        cen_evcent = centralization.evcent(go)$centralization # 特征向量中心性
      )
    }

    n_index <- apply(n_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(graph_attr(go) %>% unlist() %>% is.null())) n_index <- data.frame(graph_attr(go), n_index)
  }
  if ("v" %in% mode) {
    # Calculate Vertices Parameters
    v_index <- data.frame(
      degree = degree(go),
      clusteringC = transitivity(go, type = "local"), # 局部聚类系数
      betweenness = betweenness(go), # 节点介数
      eccentricity = eccentricity(go),
      closeness=closeness(go),
      hub_score=hub_score(go)[["vector"]]
      #page_rank = page.rank(go)$vector
      #igraph::evcent(g)[["vector"]]
      #igraph::local_efficiency(go)
    )

  if(!is.null(E(go)$cor)){
    igraph::as_data_frame(go)->edge_list
    edge_list%>%select(from,cor)%>%rbind(.,select(edge_list,to,cor)%>%rename(from=to))%>%
      group_by(from)%>%summarise(w_degree=sum(cor))->w_degree
    v_index$w_degree=w_degree[match(rownames(v_index),w_degree$from),"w_degree"]%>%unlist()
  }

    v_index <- apply(v_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(vertex_attr(go) %>% unlist() %>% is.null())) v_index <- data.frame(vertex_attr(go), v_index)
  }
  if ("e" %in% mode) {
    # Calculate Edges Parameters
    e_index <- data.frame(
      igraph::as_data_frame(go)
    )
    # if(!(edge_attr(go)%>%unlist()%>%is.null()))e_index=data.frame(edge_attr(go),e_index)

  }
  return(list(n_index = n_index, v_index = v_index, e_index = e_index))
}

#' Fit power-law distribution for a igraph
#'
#' @param go igraph
#'
#' @return ggplot
#' @export
#'
#' @examples
#' fit_power(co_net)
fit_power<-function(go){
  #度分布统计
  degree_dist <- table(degree(go))
  dat <- data.frame(degree = as.numeric(names(degree_dist)), count = as.numeric(degree_dist))
  #拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
  mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
  summary(mod)
  #提取关键值
  a <- round(coef(mod)[1], 3)
  b <- round(coef(mod)[2], 3)
  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, 3)
  #p 值可以根据置换检验的原理获得
  #将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
  #比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
  p_num <- 1
  dat_rand <- dat
  for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (999+1)

  p <- ggplot(dat, aes(x = degree, y = count)) +
    geom_point() +theme_bw() +
    stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
    labs(x = 'Degree', y = 'Count')

  #添加公式拟合的注释
  label <- data.frame(
    x=0.8*max(dat$degree),
    y=c(0.9,0.85,0.8)*max(dat$count),
    formula = c(sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
                sprintf('italic(R^2) == %.3f', R2),
                sprintf('italic(P) < %.3f', p_value))
  )

  p + geom_text(aes(x=x,y=y,label = formula), data = label,parse = T)
}

###common_characteristic
#A network can be said "smallworld" if its smallworldness is higher than one (a stricter rule is smallworldness>=3; Humphries & Gurney, 2008)
# qgraph::smallworldness(g,B = 5)
# qgraph::smallworldIndex(g)
if(F){
  # degree_distribution,KS.p>0.05 indicated fit power-law distribution.
  fit_power_law(degree(go)+1,10)
  rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  fit_power_law(degree(rand.g)+1,10)
  fit_power(go)

  #smallworldness, the average path lengths which were close to logarithms of the total number of network nodes
  #生成具有 Nv=25 个节点，r=5 个邻居，重连概率 p=0.05 的网络
  g_rand_s <- watts.strogatz.game(dim = 1, size = 25, nei = 5, p = 0.05)
  transitivity(g_rand_s)
  average.path.length(g_rand_s)
  log10(length(V(g_rand_s)))
  #具有相同节点和边数量的经典随机图
  g_rand_ER <- erdos.renyi.game(n = vcount(g_rand_s), p = ecount(g_rand_s), type = 'gnm')
  transitivity(g_rand_ER)
  average.path.length(g_rand_ER)

  #modularity,values >0.4 suggest that the network has a modular structure; Newman, 2006
  modu_dect(go)%>%graph.attributes()

  #hierarchy, R2 values of the linear relationship between logarithms of clustering coefficients and the logarithms of connectivity
  data.frame(
  y=transitivity(go,"local"),
  x=log(degree(go)))->cc_k
  summary(lm(y~x,cc_k))

}


#' Degree distribution comparison with random network
#'
#' @param go igraph object
#'
#' @return ggplot
#' @export
#'
#' @examples
#'
#' rand_net(co_net)
#' rand_net_par(co_net)->a
rand_net<-function(go = go){
  lib_ps("igraph")
  #generate a random network
  rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")

  data1 = data.frame(freq= degree_distribution(go),net = "network",degree = 1:length(degree_distribution(go)))
  data2 = data.frame(freq = degree_distribution(rand.g) ,net = "random E-R",degree = 1:length(degree_distribution(rand.g)))
  data = rbind(data1,data2)
  p1 <- ggplot(data)+
    geom_point(aes(x = degree,y = freq,group =net,fill = net),pch = 21,size = 2) +
    geom_smooth(aes(x = degree,y = freq,group =net,color = net),se=F)+
    labs(x="Degree",y="Proportion")+scale_color_manual(values = c("#F58B8B","#7AADF0"))+
    scale_fill_manual(values = c("#F58B8B","#7AADF0"))+
    theme_pubr()+theme(legend.position = c(0.8,0.9),legend.title = element_blank())
  return(p1)
}

#' Net_pars of many random network
#'
#' @param go igraph
#' @param reps simulation time
#' @param threads threads
#'
#' @export
rand_net_par<-function(go,reps=99,threads=4){
  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  rand_net_pars <- foreach(i = 1:reps,.options.snow = opts,
                          .packages = c("igraph"),.combine = "rbind") %dopar% {
                            #generate a random network
                            rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
                            pctax::net_par(rand.g,mode = "n")[["n_index"]]->indexs
                            wc <- igraph::cluster_fast_greedy(rand.g)
                            indexs$modularity <- modularity(wc)
                            indexs
                          }
  stopCluster(cl)
  if(F){
    ggplot(a,aes(x=clusteringC))+geom_histogram()+
      geom_vline(xintercept = transitivity(go), col = 'red')
    ggplot(a,aes(x=ave_path_len))+geom_histogram()+
      geom_vline(xintercept = average.path.length(go), col = 'red')+xlim(2,3)
  }
  rand_net_pars
}

#' Calculate small-world coefficient
#'
#' @param go igraph
#'
#' @export
#'
smallworldness<-function(go){
  rand_net_par(go,reps = 99)->rands
  small_world_coefficient=(transitivity(go)/mean(rands$clusteringC))/
    (average.path.length(go)/mean(rands$ave_path_len))
  small_world_coefficient
}

#========6.modules=========
#' Detect the modules
#'
#' @param go a igraph object
#' @param method cluster_method:("cluster_walktrap","cluster_edge_betweenness","cluster_fast_greedy","cluster_spinglass")
#' @return a igraph object
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

  graph_attr(go)$modularity<- modularity(wc)

  rand.g <- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  rand_m=modularity(cluster_fast_greedy(rand.g))
  r_modularity=(modularity(wc)-rand_m)/rand_m#相对模块指数
  graph_attr(go)$r_modularity=r_modularity
  return(go = go)
}

#' Zi-Pi calculate
#'
#' @param go igraph object after modu_dect()
#' @param mode use 7-group (mode=1) or 4-group (mode=2)
#'
#' @return igraph
#' @export
#'
#' @references 1. Guimerà, R. & Amaral, L. Functional cartography of complex metabolic networks. (2005) doi:10.1038/nature03288.
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu,mode = 2)->co_net_modu
#' zp_plot(co_net_modu)
zp_analyse<-function(go,mode=1){
  # =====模块内模块间分析
  v_index=igraph::as_data_frame(go,"vertices")
  if(!"module"%in%names(v_index))stop("no modules, please modu_dect() first")
  if("roles"%in%names(v_index))message("areadly has roles, overwrite!")
  within <- within_module_deg_z_score(go)
  v_index$Ki=within$Ki;v_index$Zi=within$Zi
  pc <- part_coeff(go)
  v_index$Pi=pc$Pi

  if(mode==1){
    lab <- c("Ultra-peripherals","Peripherals","Non-hub connectors","Non-hub kinless nodes","Provincial hubs", "Connector hubs","Kinless hubs")
    backs=data.frame(
      x1= c(0,0.05, 0.62,0.8, 0,0.3,0.75),
      x2= c(0.05,0.62,0.8, 1, 0.3,0.75, 1),
      y1= c(-Inf,-Inf,-Inf,-Inf, 2.5, 2.5,2.5),
      y2= c(2.5, 2.5, 2.5,2.5, Inf, Inf, Inf),
      lab=factor(lab,levels = lab)
    )
  }
  else if(mode==2){
    lab <- c("Peripherals", "Network hubs", "Module hubs","Connectors")
    backs=data.frame(
      x1= c(0, 0.62, 0, 0.62),
      x2= c(0.62, 1, 0.62, 1),
      y1= c(-Inf, 2.5, 2.5, -Inf),
      y2= c(2.5, Inf, Inf, 2.5),
      lab=factor(lab,levels = lab)
    )
  }
  deter_role=\(x,y,backs=backs){
    for (i in 1:nrow(backs)){
      if(between(x,backs$x1[i],backs$x2[i])&&between(y,backs$y1[i],backs$y2[i])){
        role=backs$lab[i]
        break
      }
    }
    return(role)
  }
  v_index$roles=apply(v_index, 1, \(x)deter_role(x["Pi"],x["Zi"],backs))
  vertex.attributes(go)<-as.list(v_index)
  return(go)
}

#' Plot a community (modularity)
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
#' modu_plot(co_net_modu, n_modu = 50)
modu_plot <- function(go, coors=NULL,n_modu = 1, ...) {
  if(!"module"%in%vertex_attr_names(go))stop("no modules, please modu_dect() first")
  as_group <- \(x){
    x <- list(membership = x)
    vids <- names(x$membership)
    modus <- tapply(vids, x$membership, simplify = FALSE, function(x) x)
    return(modus)
  }

  set.seed(123)
  if (is.null(coors)) coors <- layout.fruchterman.reingold(go)
  else {
    coors<-coors[match(V(go)$name,coors$name),2:3]%>%as.matrix()
  }
  members=V(go)$module;names(members)=V(go)$name

  if (n_modu) {
  # 筛选成员大于n_modu的模块
  (which((members %>% table()) > n_modu))->cho_modu
  members[members %in%cho_modu ] %>% as_group() -> new_modu
  col=ifelse(members %in%cho_modu,members ,"grey")

  plot(go,...,
    main = "Modularity network", layout = coors,
    mark.groups = new_modu, vertex.color = col,
    vertex.label.font = 2, vertex.label.color = "black",
    vertex.label.cex = .05* V(go)$size, vertex.size = V(go)$size,
    vertex.label = ifelse(V(go)$size > 8, V(go)$name, NA),
    edge.lty = 1, edge.curved = TRUE, margin = c(0, 0, 0, 0)
  )
  leg1<-paste(names(new_modu),sapply(new_modu, length),"nodes",sep = " :")
  legend(1, 1, cex = .7, legend =leg1 ,
         fill = rainbow(length(new_modu),alpha = 0.3), bty = "n", title = "Modules")

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
#' @param go igraph object after zp_analyse()
#' @param label show label or not
#'
#' @return a ggplot object
#' @export
#' @import ggpubr ggrepel
#' @seealso \link{zp_analyse}
zp_plot <- function(go,label=T) {
  lib_ps("ggrepel","ggpubr","igraph")
  as.data.frame(vertex.attributes(go))->taxa.roles
  if(!"roles"%in%names(taxa.roles))stop("no roles, please zp_analyse() first")
  mode=ifelse(nlevels(V(go)$roles)==7,1,2)
  if(mode==1){
    lab <- c("Ultra-peripherals","Peripherals","Non-hub connectors","Non-hub kinless nodes","Provincial hubs", "Connector hubs","Kinless hubs")
    CPCOLS<-c("#FCF6EFFC","#EEBCF5", "#EDEDA4", "#FAA371","#FC5D6096" ,"#9BC799B9" ,"#94CCF2AC")
    names(CPCOLS)<-lab
    backs=data.frame(
      x1= c(0,0.05, 0.62,0.8, 0,0.3,0.75),
      x2= c(0.05,0.62,0.8, 1, 0.3,0.75, 1),
      y1= c(-Inf,-Inf,-Inf,-Inf, 2.5, 2.5,2.5),
      y2= c(2.5, 2.5, 2.5,2.5, Inf, Inf, Inf),
      lab=factor(lab,levels = lab)
    )
  }
  else if(mode==2){
    lab <- c("Peripherals", "Network hubs", "Module hubs","Connectors")
    CPCOLS<-c("#FCF6EFFC" ,"#FC5D6096" ,"#9BC799B9" ,"#94CCF2AC")
    names(CPCOLS)<-lab
    backs=data.frame(
      x1= c(0, 0.62, 0, 0.62),
      x2= c(0.62, 1, 0.62, 1),
      y1= c(-Inf, 2.5, 2.5, -Inf),
      y2= c(2.5, Inf, Inf, 2.5),
      lab=factor(lab,levels = lab)
    )
  }

  p <- ggplot() +
    geom_rect(data = backs, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab)) +
    guides(fill = guide_legend(title = "Topological roles")) +
    scale_fill_manual(values = CPCOLS) +
    geom_point(data = taxa.roles, aes(x = Pi, y = Zi, color = factor(class))) +
    ggpubr::theme_pubr(legend = "right") +
    guides(colour = "none") +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Participation Coefficient") +
    ylab("Within-module connectivity z-score")

  if(label){
    p=p+ggrepel::geom_text_repel(
        data = filter(taxa.roles, !roles%in%c("Peripherals","Ultra-peripherals")),
        aes(x = Pi, y = Zi, label = name)
      )
  }
  return(p)
}

#========7.stability===========
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
#' robustness_test(co_net,step=4)->robust_res
#' class(robust_res)
#' plot(robust_res,index="ave_degree",mode=2)
robustness_test <- function(go, partial = 0.5,step=10,reps=9,threads=4) {
  lib_ps("igraph")
  cal_del<-\(go,partial,step,rep){
    nodes <- length(V(go))
    floor(nodes * partial) -> del_i
    del_i_indexs <- data.frame()
    sequ=seq(0,del_i,step)
    if(sequ[length(sequ)]<del_i)sequ=c(sequ,del_i)

    for (i in sequ) {
      # 在网络中随机移除 i 个节点
      remove_node <- sample(1:nodes, i)
      dp <- igraph::delete.vertices(go, remove_node)
      dp=delete.vertices(dp,V(dp)[degree(dp)==0])

      # 计算指标
      tmp_ind=(pctax::net_par(dp,mode = "n")$n_index)
      del_i_indexs <- rbind(
        del_i_indexs,
        data.frame(tmp_ind, i = i)
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
                                         .packages = c("igraph")) %dopar%{
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
#' @param mode plot mode
#' @param indexs indexs selected to show
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return a ggplot
#' @method plot robustness
#' @export
#' @rdname robustness_test
plot.robustness<-function(robust_res, indexs = c("nat_connectivity", "ave_path_len", "ave_degree"),mode=1){
  lib_ps("reshape2")
  robust_res %>%
    dplyr::select(i, all_of(indexs)) %>%
    reshape2::melt(id.var = "i") -> pdat

  pdat%>%group_by(i,variable)%>%summarise(mean=mean(value),sd=sd(value),se=sd/sqrt(length(value)))->sdd
  if(mode==1){
  p <- ggplot(sdd, aes(x=i, y=mean,col=indexs)) +
    geom_line()+
    geom_errorbar(data = sdd,aes(ymax=mean+se,ymin=mean-se))+
    #geom_smooth(se = FALSE,method = "loess",formula = 'y ~ x') +
    facet_wrap(. ~ variable, scales = "free")
  }
  if(mode==2){
  p<-ggplot(sdd,aes(x=i, y=mean,col=indexs))+
    geom_point(size=0.2,alpha=0.4)+
    geom_smooth(se = F,method = "loess",formula = 'y ~ x')+
    ggpmisc::stat_poly_eq(
      aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = '~~~~~')),
      formula = y ~ x,  parse = TRUE,label.x = "right",
      size = 3, #公式字体大小
    )+
    facet_wrap(. ~ variable, scales = "free")
  }

  return(p)
}

#https://mp.weixin.qq.com/s?__biz=MzAwODk1Njk5MA==&mid=2247485401&idx=1&sn=30f4c324e1013482b98aa5b5b4f5cdc1&scene=21#wechat_redirect

#' Cohesion calculate
#'
#' @param otutab otutab
#' @param reps iteration time
#' @param threads threads
#' @param mycor a correlation matrix you want to use, skip the null model build
#'
#' @return a list with two dataframe
#' @export
#' @references 1. Herren, C. M. & McMahon, K. Cohesion: a method for quantifying the connectivity of microbial communities. (2017) doi:10.1038/ismej.2017.91.
#' @examples
#' Cohesion(otutab[1:50,])->a
#' stackplot(abs(t(a$Cohesion)),metadata,groupID = "Group")
Cohesion<-function(otutab,reps= 200,threads=4,mycor=NULL){
  d <- t(otutab)
  rel.d <- d / rowSums(d)

  if(is.null(mycor)){# Create observed correlation matrix
  cor.mat.true <- cor(rel.d)
  lib_ps("foreach")
  lib_ps("doSNOW")
  nc <- ncol(rel.d)
  pb <- txtProgressBar(max =nc, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  a=foreach (which.taxon = 1:nc,.options.snow=opts) %dopar%{
    #create vector to hold correlations from every permutation for each single otu
    ## perm.cor.vec.mat stands for permuted correlations vector matrix
    perm.cor.vec.mat <- vector()

    for(i in 1:reps){
      #Create empty matrix of same dimension as rel.d
      perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
      rownames(perm.rel.d) <- rownames(rel.d)
      colnames(perm.rel.d) <- colnames(rel.d)

      #For each otu
      for(j in 1:dim(rel.d)[2]){
        # Replace the original taxon vector with a permuted taxon vector
        perm.rel.d[, j ] <- sample(rel.d[ ,j ])
      }

      # Do not randomize focal column
      perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]

      # Calculate correlation matrix of permuted matrix
      cor.mat.null <- cor(perm.rel.d,perm.rel.d[, which.taxon])

      # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null)

    }

    # Save the median correlations between the focal taxon and all other taxa
    apply(perm.cor.vec.mat, 1, median)
  }
  simplify2array(a)->med.tax.cors
  stopCluster(cl)

  obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
  else obs.exp.cors.mat=mycor

  diag(obs.exp.cors.mat) <- 0
  # Calculate connectedness by averaging positive and negative observed - expected correlations
  pn.mean=function(vector,mode="p"){
    if(mode=="p")pos.vals=vector[which(vector > 0)]
    else pos.vals=vector[which(vector < 0)]
    p.mean <- mean(pos.vals)
    if(length(pos.vals) == 0) p.mean <- 0
    return(p.mean)
  }
  connectedness.pos <- apply(obs.exp.cors.mat, 2, pn.mean,mode="p")
  connectedness.neg <- apply(obs.exp.cors.mat, 2, pn.mean,mode="n")

  # Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
  cohesion.pos <- rel.d %*% connectedness.pos
  cohesion.neg <- rel.d %*% connectedness.neg

  ####
  #### Combine vectors into one list and print
  output <- list(data.frame(neg=connectedness.neg, pos=connectedness.pos),
                 data.frame(neg=cohesion.neg, pos=cohesion.pos))
  names(output) <- c("Connectedness", "Cohesion")
  return(output)
}

#' Vulnerability
#'
#' @param go igraph
#'
#' @return a vector
#' @export
#' @description
#' \deqn{Vi=\frac{E-Ei}{E}}
#' E is the global efficiency and Ei is the global efficiency after the removal of the node i and its entire links.
vulnerability<-function(go,threads=1){
  if(is.null(V(go)$name))V(go)$name=V(go)
  if(threads==1){
    lapply(V(go)$name,\(i)global_efficiency(delete_vertices(go,i)))->tmpls
    do.call(c,tmpls)->tmpv
  }
  else{
    #parallel
    lib_ps("foreach","doSNOW")
    nc=length(V(go)$name)
    pb <- txtProgressBar(max =nc, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    cl <- makeCluster(threads)
    registerDoSNOW(cl)
    tmpv <- foreach(i = V(go)$name,.combine = "c",
                    .options.snow = opts,
                    .packages = c("igraph"))%dopar%{
                      global_efficiency(delete_vertices(go,i))
                    }
    stopCluster(cl)
  }
  (global_efficiency(go)-tmpv)/global_efficiency(go)
}

#' Robustness after remove 50\% nodes or some hubs
#'
#' @param go igraph
#' @param keystone remove 70\% keystone (default:False)
#' @param reps simulation time
#' @param threads threads
#'
#' @export
#'
#' @examples
#' robustness(co_net)
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu,mode = 2)->co_net_modu
#' robustness(co_net_modu,keystone=T)
robustness<-function(go,keystone=F,reps=9,threads=4){
  nodes <- length(V(go))
  floor(nodes * 0.5) -> del_i
  if(keystone){
    as.data.frame(vertex.attributes(co_net_modu))->tmp_v
    tmp_v%>%filter(roles=="Module hubs")%>%pull(name)->hubs
    floor(length(hubs)* 0.7) -> del_i
  }

  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  del_i_indexs<- foreach (i = 1:reps,.combine = "rbind",
                          .options.snow = opts,
                          .packages = c("igraph","dplyr")) %dopar%{
                            if(!keystone)remove_node <- sample(1:nodes, del_i)
                            if(keystone)remove_node <- sample(hubs, del_i)

                            dp <- igraph::delete.vertices(go, remove_node)
                            dead_s="init"
                            #calculated the abundance-weighted mean interaction strength (wMIS) of nodes,<=0 will dead
                            #repeat until all >0
                            while(length(dead_s)>0){
                              igraph::as_data_frame(dp)->edge_list
                              edge_list%>%select(from,cor)%>%rbind(.,select(edge_list,to,cor)%>%rename(from=to))%>%
                                group_by(from)%>%summarise(w_degree=sum(cor))->w_degree
                              w_degree%>%filter(w_degree<=0)%>%pull(from)->dead_s
                              dp=delete.vertices(dp,dead_s)
                            }

                            # 计算指标
                            tmp_ind=(pctax::net_par(dp,mode = "n")$n_index)
                            data.frame(tmp_ind, reps = i)
                          }
  stopCluster(cl)
  return(del_i_indexs)
}

