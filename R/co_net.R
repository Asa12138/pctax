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
#' sample(t(otutab) %>% as.data.frame(), 50) -> totu2
#' colnames(totu2) <- paste0("tt", 1:50)
#' c_net_cal(totu, totu2) -> corr2
c_net_cal <- function(totu, totu2 = NULL, filename = "occor",threads=4) {
  corr<-par_cor(totu,totu2,threads = threads,method = "spearman",p.adjust.methods=NULL)
  r <- round(corr$r, 5)
  p.value <- round(corr$p.value, 5)
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
  if(is.null(totu2)){
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=nc)
      for(j in 1:nc) {
        if(j > i) corr1[,j] <- r(totu[,i],totu[,j])
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
        corr1[,j] <- r(totu[,i],totu2[,j])
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
  if(!is.null(edge.attributes(up)$weight))up<-delete_edge_attr(up,"weight")
  if ("n" %in% mode) {
    # Calculate Network Parameters
    num_nodes <- length(V(go)) # number of nodes
    num_edges <- length(E(go)) # number od edges
    edge_density <- edge_density(go) # density of network, connectance
    ave_path_len <- average.path.length(up) # Average path length
    #w_ave_path_len <- ifelse(is.null(E(go)$weight), ave_path_len, average.path.length(go)) # weighted Average path length
    ave_degree <- mean(igraph::degree(go)) # Average degree
    w_ave_degree <- ifelse(is.null(E(go)$weight), ave_degree, sum(E(go)$weight) / length(V(go))) # weighted degree
    diameter <- diameter(up) # network diameter
    #w_diameter <- ifelse(is.null(E(go)$weight), diameter, diameter(go))
    clusteringC <- transitivity(go) # Clustering coefficient
    cen_betweenness <- centralization.betweenness(go)$centralization # 介数中心性(Betweenness centralization)
    nat_connectivity <- nc(go) # 自然连通度
    neg_percent<-sum(E(go)$cor<0)/num_edges # 负相关比例

    n_index <- data.frame(num_nodes, num_edges, edge_density, neg_percent,ave_path_len, ave_degree,
       w_ave_degree,diameter,clusteringC,cen_betweenness,nat_connectivity)
    if(!fast){
      # mean_dist=mean_distance(go)#平均距离
      # w_mean_dist=ifelse(is.null(E(go)$weight),mean_dist,mean_distance(go))
      # v_conn <- vertex.connectivity(go) # 点连通度
      # e_conn <- edge.connectivity(go) # 边连通度
      # components <- count_components(go) # 连通分量的数目
      cen_closeness <- centralization.closeness(go)$centralization # 紧密中心性
      cen_degree <- centralization.degree(go)$centralization # 度中心性(Degree centralization)
      cen_evcent <- centralization.evcent(go)$centralization # 特征向量中心性

      n_index <- data.frame(
        num_nodes, num_edges, edge_density,neg_percent, ave_path_len, ave_degree, w_ave_degree,
        diameter, clusteringC, cen_betweenness,nat_connectivity, cen_closeness, cen_degree, cen_evcent
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
      page_rank = page.rank(go)$vector
    )
  if(!is.null(E(go)$weight)){
    igraph::as_data_frame(go)->edge_list
    edge_list%>%select(from,weight)%>%rbind(.,select(edge_list,to,weight)%>%rename(from=to))%>%
      group_by(from)%>%summarise(w_degree=sum(weight))->w_degree
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
#' modu_plot(co_net_modu, n_modu = 5)
modu_plot <- function(wc, go = NULL, coors=NULL,n_modu = 1, ...) {
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
  else {
    coors<-coors[match(V(go)$name,coors$name),2:3]%>%as.matrix()
  }
  if (n_modu) {
  # 筛选成员大于n_modu的模块
  (which((membership(wc) %>% table()) > n_modu))->cho_modu
  membership(wc)[membership(wc) %in%cho_modu ] %>% as_group() -> new_modu
  col=ifelse(membership(wc) %in%cho_modu,membership(wc) ,"grey")

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
  lib_ps("ggrepel","ggpubr")
  taxa.roles <- indexs$v_index
  x1 <- c(0, 0.62, 0, 0.62)
  x2 <- c(0.62, 1, 0.62, 1)
  y1 <- c(-Inf, 2.5, 2.5, -Inf)
  y2 <- c(2.5, Inf, Inf, 2.5)
  lab <- c("Peripheral", "Network hubs", "Module hubs","Connectors")
  lab=factor(lab,levels = lab)
  CPCOLS <- c( "#FCF6EFFC","#FC5D6096", "#9BC799B9", "#94CCF2AC")
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
#' robustness_test(co_net,step=4)->robust_res
#' class(robust_res)
#' plot(robust_res,index="degree",mode=2)

robustness_test <- function(go, partial = 0.7,step=1,reps=9,threads=4) {
  lib_ps("igraph")
  cal_del<-\(go,partial,step,rep){
    nodes <- length(V(go))
    floor(nodes * partial) -> del_i
    del_i_indexs <- data.frame()

    for (i in seq(0,del_i,step)) {
      # 在网络中随机移除 i 个节点
      remove_node <- sample(1:nodes, i)
      dp <- igraph::delete.vertices(go, remove_node)
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

#========Cohesion======
#https://mp.weixin.qq.com/s?__biz=MzAwODk1Njk5MA==&mid=2247485401&idx=1&sn=30f4c324e1013482b98aa5b5b4f5cdc1&scene=21#wechat_redirect

