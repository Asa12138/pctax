#phylogenic tree==============
#
taxaclass=c( "Kingdom","Phylum","Class","Order","Family","Genus","Speies" ,"Rank8", "Rank9" ,"Rank10")

#' Complete a taxonomy table
#'
#' @param taxdf taxonomy table
#' @param type species
#' @import dplyr tibble
#' @return a good taxonomy table
#' @export
#' @references MicrobiotaProcess
#' @examples
#'taxmat = matrix(sample("onelevel", 7*2, replace = TRUE), nrow = 2, ncol = 7)%>%as.data.frame()
#'colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#'fillNAtax(taxmat)
fillNAtax<-function (taxdf, type = "species") {
  lib_ps("dplyr","tibble")

  grepl.data.frame<-\(pattern, x, ...) {
    y <- if (length(x)) {
      do.call("cbind", lapply(x, "grepl", pattern = pattern,...))
    }
    else {
      matrix(FALSE, length(row.names(x)), 0)
    }
    if (.row_names_info(x) > 0L)
      rownames(y) <- row.names(x)
    y
  }

  remove_unclassfied<-\(taxdf) {
    taxdf[grepl.data.frame("Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig",
                           taxdf, ignore.case = TRUE)] <- NA
    return(taxdf)
  }

  remove_na_taxonomy_rank<-function (x) {
    x <- as.data.frame(x, check.names = FALSE)
    x <- remove_unclassfied(x)
    indx <- vapply(x, function(i) all(is.na(i)), FUN.VALUE = logical(1)) %>%
      setNames(NULL)
    x <- x[, !indx, drop = FALSE]
    return(x)
  }

  taxlevelchar<-c("k","p","c","o","f","g","s","st")

  addtaxlevel<-function (taxdf) {
    taxlevelchar <- taxlevelchar[seq_len(length(taxdf))]
    paste(taxlevelchar, taxdf, sep = "__")
  }

  duplicatedtaxcheck<-function (taxdf) {
    if (ncol(taxdf) == 1) {
      return(taxdf)
    }
    taxdf <- taxdf %>% rownames_to_column()
    for (i in ncol(taxdf):3) {
      tmp <- split(taxdf, taxdf[, i])
      for (j in seq_len(length(tmp))) {
        flag <- length(unique(as.vector(tmp[[j]][, i - 1])))
        if (flag > 1) {
          tmp[[j]][, i] <- paste(tmp[[j]][, i], tmp[[j]][,
                                                         i - 1], sep = "_")
        }
      }
      taxdf <- do.call("rbind", c(tmp, make.row.names = FALSE))
    }
    return(taxdf)
  }

  repduplicatedtaxcheck<-function (taxdf) {
    for (i in seq_len(7)) {
      taxdf <- duplicatedtaxcheck(taxdf) %>% column_to_rownames(var = "rowname")
    }
    return(taxdf)
  }

  filltaxname<-function (taxdf, type = "species") {
    tmprownames <- rownames(taxdf)
    #about 2chars more
    indexmark <- apply(taxdf, 2, function(x) {
      nchar(x, keepNA = TRUE)
    }) <= 4
    taxdf[indexmark] <- NA
    if (any(is.na(taxdf[, 1]))) {
      if (type == "species") {
        prefix <- "k__"
      }
      else {
        prefix <- "d1__"
      }
      taxdf[is.na(taxdf[, 1]), 1] <- paste0(prefix, "Unknown")
    }
    indextmp <- apply(is.na(taxdf), 1, which)
    if (length(indextmp) == 0) {
      taxdf <- data.frame(taxdf, check.names = FALSE)
      return(taxdf)
    }
    taxdf <- apply(taxdf, 1, zoo::na.locf)
    taxdf <- lapply(seq_len(ncol(taxdf)), function(i) taxdf[,
                                                            i])
    taxdf <- data.frame(t(mapply(newtaxname, taxdf, indextmp)),
                        stringsAsFactors = FALSE)
    rownames(taxdf) <- tmprownames
    return(taxdf)
  }

  newtaxname<-function (x, y) {
    y <- as.vector(y)
    x[y] <- paste(taxlevelchar[y], x[y], sep = "__un_")
    x
  }

  rown=rownames(taxdf)
  taxdf <- remove_na_taxonomy_rank(taxdf)
  if (type != "species") {
    assign("taxlevelchar", paste0("d", seq_len(ncol(taxdf))),
           envir = .GlobalEnv)
  }
  else {
    assign("taxlevelchar", c("k", "p", "c", "o", "f", "g",
                             "s", "st"), envir = .GlobalEnv)
  }
  if (!(grepl("^k__", taxdf[1, 1]) || grepl("^d1__", taxdf[1, 1]))) {
    tmprownames <- rownames(taxdf)
    tmpcolnames <- colnames(taxdf)
    taxdf <- t(apply(taxdf, 1, as.character))
    taxdf[is.na(taxdf)] <- ""
    taxdf <- data.frame(t(apply(taxdf, 1, addtaxlevel)),
                        stringsAsFactors = FALSE)
    rownames(taxdf) <- tmprownames
    colnames(taxdf) <- tmpcolnames
  }
  taxdf <- filltaxname(taxdf, type = type)
  taxdf <- repduplicatedtaxcheck(taxdf)
  taxdf<-taxdf[rown,]
  attr(taxdf, "fillNAtax") <- TRUE
  #taxdf%>%mutate_all(.funs = \(x)gsub(" ","_",x))->taxdf
  return(taxdf)
}



#' Before df2tree check
#'
#' @param f_tax table
#'
#' @return table
#' @export
#'
before_tree<-function(f_tax){
  for (i in seq_len(ncol(f_tax)-1)) {
    tmp=distinct(f_tax[,c(i,i+1)])
    du=tmp[duplicated(tmp[,2]),2]
    if(length(du)>0){
      print(du)
      f_tax[,i+1]=apply(f_tax, 1, \(x)ifelse(x[i+1]%in%du,paste(x[i+1],x[i],sep = "_"),x[i+1]))
    }
  }
  f_tax
}
#' From a dataframe to construct a phylo
#'
#' @param taxa dataframe
#'
#' @return phylo object
#' @export
#'
#' @examples
#' data(otutab)
#' df2tree(taxonomy)->tax_tree
#' #check all nodes matched!
#' picante::match.phylo.comm(tax_tree,t(otutab))->nn
#' nrow(nn$comm)==nrow(t(otutab))
#'
df2tree<-function(data){
  data <- data.frame(Root = rep("r__root", nrow(data)), data)
  datalist <- list()
  clnm <- colnames(data)
  for (i in seq_len(ncol(data) - 1)) {
    tmpdat <- data[, c(i, i + 1)]
    colnames(tmpdat) <- c("parent", "child")
    tmpdat %<>% dplyr::mutate(nodeClass = clnm[i + 1], nodeDepth = i) %>%
      dplyr::distinct()
    datalist[[i]] <- tmpdat
  }
  datalist <- do.call("rbind", datalist)
  datalist <- datalist[!duplicated(datalist), ]
  isTip <- !as.vector(datalist$child) %in% as.vector(datalist$parent)
  index <- rep(NA, length(isTip))
  index[isTip] <- seq(1, sum(isTip))
  index[!isTip] <- seq(sum(isTip) + 2, length(isTip) + 1)
  mapping <- data.frame(node = index, labelnames = as.vector(datalist$child),
                        isTip)
  indxx <- match(mapping$labelnames, datalist$child)
  mapping$nodeClass <- datalist[indxx, "nodeClass"]
  mapping$nodeDepth <- datalist[indxx, "nodeDepth"]
  parentnode <- mapping[match(as.vector(datalist$parent),
                              as.vector(mapping$labelnames)), ]$node
  childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)),
  ]$node
  edges <- cbind(parentnode, childnode)
  colnames(edges) <- NULL
  edges[is.na(edges)] <- sum(isTip) + 1
  root <- data.frame(node = sum(isTip) + 1, labelnames = "r__root",
                     isTip = FALSE, nodeClass = "Root", nodeDepth = 0)
  mapping <- rbind(root, mapping)
  mapping <- mapping[order(mapping$node), ]
  node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
  tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
  mapping <- mapping[, colnames(mapping) %in% c("node", "nodeClass",
                                                "nodeDepth")]
  taxphylo <- structure(list(edge = edges, node.label = node.label,edge.length=rep(1,nrow(edges)),
                             tip.label = tip.label, Nnode = length(node.label)),
                        class = "phylo")
  return(taxphylo)
}

#会将带空格的label改变，淘汰
makeNewick<-function (taxa) {
  lib_ps("ggtree")
  #taxa%>%mutate_all(.funs = \(x)gsub(" ","",x))->taxa
  makeNewick1<-\(taxa, naSub = "_") {
    if (!is.null(naSub))
      taxa[is.na(taxa)] <- naSub
    if (ncol(taxa) == 0)
      return("")
    bases <- unique(taxa[, 1])
    innerTree <- sapply(bases, function(ii) makeNewick1(taxa[taxa[, 1] == ii, -1, drop = FALSE]))
    out <- sprintf("(%s)", paste(sprintf("%s%s", innerTree,bases), collapse = ","))
    return(out)
  }
  paste0(makeNewick1(taxa),';')->tree
  tree%>%read.tree(text = .)->nwk
  nwk$edge.length=rep(1,nrow(nwk$edge))
  return(nwk)
}

#' Annotate a tree
#'
#' @param f_tax taxonomy dataframe
#' @param otutab otutab, rowname==rowname(taxonomy)
#' @param level 1~7
#'
#' @return a treedata
#' @export
#'
#' @examples
#' data(otutab)
#' ann_tree(taxonomy,otutab)->tree
#' easy_tree(tree)
ann_tree<-function(f_tax,otutab,level=ncol(f_tax)){
  if(any(rownames(f_tax)!=rownames(otutab)))stop("rowname not match")
  otutab%>%rowSums()->num
  res=data.frame(label="",abundant=sum(num))
  for (i in 1:level){
    aggregate(num,by=list(f_tax[,i]),sum)->tmp
    colnames(tmp)<-colnames(res)
    res=rbind(res,tmp)
  }
  #le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  le=c("root",colnames(f_tax))
  df2tree(f_tax[,1:level])%>%ggtree::fortify()->tree
  tree$level<-le[ceiling(tree$branch)+1]
  #tree$level<-factor(tree$level,levels = le)
  left_join(tree,tree[,c("node","label")],by=c("parent"="node"))%>%dplyr::rename(label=label.x,parent_label=label.y)->tree1
  left_join(tree1,res)->tree2
  ann_tax=data.frame()
  for(i in level:1){
    f_tax[,1:i,drop=F]%>%distinct()->tmpdf
    tmpdf[,ncol(tmpdf)]->rownames(tmpdf)
    ann_tax=vctrs::vec_c(ann_tax,tmpdf)
  }
  left_join(tree2,ann_tax%>%tibble::rownames_to_column("label"),by="label")->tree3
  return(tree3)
}

#' Easy way to plot a phylogenetic tree
#'
#' @param tree result from ann_tree
#'
#' @return a ggplot
#' @export
#'
#' @rdname ann_tree
easy_tree<-function(tree,highlight="Phylum",colorfill="color"){
  if(!requireNamespace("treedataverse"))BiocManager::install("YuLab-SMU/treedataverse");
  library(treedataverse)
  #可以剪切部分树
  tree%>%mutate(label1=sub(".?__","",label))->tree2
  filter(tree2,level==highlight)%>%select(node,label1)->phy
  fontsize=round(600/sum(tree2$isTip),2)
  if(fontsize>5)fontsize=5
  if(fontsize<0.2)fontsize=0.2
if(colorfill=="fill"){
  ggtree(tree2,layout = 'radial',size=0.5)+
    geom_highlight(data = phy,aes(node = node, fill = label1),alpha=0.5)+
    scale_fill_manual(name=highlight,values =get_cols(nrow(phy)),na.value = NA)+
    ggnewscale::new_scale_fill()+
    #scale_fill_d3()+
    geom_fruit(
      geom=geom_tile,
      mapping=aes(y=label, fill=log(abundant)),
      stat="identity",offset = 0.05,pwidth = 0.1,
    )+scale_fill_gradient(name="log(cpm)",low = "#FFFFFF",high = "red")+
    geom_tiplab(aes(label=label1),color="black",size=fontsize,offset = 1, show.legend = FALSE)}

if(colorfill=="color"){
  tree2=groupClade(tree2,setNames(phy$node,phy$label1))
  tree2$group[tree2$group==0]=NA

  ggtree(tree2, aes(color=group),layout = 'radial',size=0.5)+
    scale_color_manual(name=highlight,values = get_cols(nrow(phy)),na.translate=FALSE)+
    guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))+
    #scale_fill_d3()+
    geom_fruit(
      geom=geom_tile,
      mapping=aes(y=label, fill=log(abundant)),
      stat="identity",offset = 0.05,pwidth = 0.1,
    )+scale_fill_gradient(name="log(cpm)",low = "#FFFFFF",high = "red")+
    geom_tiplab(aes(label=label1),color="black",size=fontsize,offset = 1, show.legend = FALSE)
  }
}

#' Get strip for a circle tree
#'
#' @param name
#' @param tree2
#' @param name_level
#' @param flat_n
#'
#' @export
#'
get_strip=\(name,tree2,flat_n=5){
  mx=max(tree2$x)
  id=filter(tree2,label==name)%>%pull(node)
  offspring(tree2,id)%>%filter(x==mx)%>%arrange(angle)->tmp
  offset.text=ifelse(nrow(tmp)>flat_n,0.5,0)
  hjust=ifelse(nrow(tmp)>flat_n,0.5,0)
  angle=ifelse(nrow(tmp)>flat_n,mean(tmp$angle)-90,mean(tmp$angle))
  list(head(tmp,1)%>%pull(label),tail(tmp,1)%>%pull(label),offset.text,hjust,angle)
}


#' Plot a sankey
#'
#' @param tree ann_tree result
#' @param top_N each level has top_N
#' @param notshow some words you don't want to show
#'
#' @export
#'
#' @examples
#' data(otutab)
#' ann_tree(taxonomy[,c(1,5,6,7)],otutab)->tree
#' sangji_plot(tree)
sangji_plot<-function(tree,top_N=5,notshow=c(),intermediate=F,width = 3000,height = 500,...){
  if(!requireNamespace("sankeyD3"))remotes::install_github("fbreitwieser/sankeyD3");
  library(sankeyD3)
  if(F){plot_ly(
    type = 'sankey', orientation = 'h',
    #节点颜色按预定义的变量的分类着色
    node = list(
      label = nodes$label,
      pad = 50, thickness = 20,
      line = list(color = 'black', width = 0.5)
    ),
    #指定预定义的变量 id，按相关性方向赋值连线颜色，相关性强度赋值连线尺寸
    link = list(
      source = links$IDsource, target =links$IDtarget,
      value = links$abundant, label = links$abundant
    )
  )}
  #桑基图
  le=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

  #select show part
  lib_ps("tidytree")
  if(length(notshow)>0)tree1=tree[!grepl(paste0(notshow,collapse = "|"),tree$label),]
  else tree1=tree
  tree1%>%group_by(level)%>%top_n(top_N,abundant)%>%ungroup()->sangji_dat

  if(intermediate){
    all_node=sapply(sangji_dat$node, \(i)tidytree::ancestor(tree,i)%>%pull(node))%>%unlist()
    all_node=c(sangji_dat$node,all_node)%>%unique()
    tree%>%filter(node%in%all_node)->sangji_dat
  }

  #show levels
  sangji_dat%>%filter(!label%in%c("","r__root"))%>%select(label,level,abundant)%>%as.data.frame()->nodes

  #tree structure
  sangji_dat%>%filter(parent_label!="r__root")%>%select(label,parent_label,abundant)->links
  links=as.data.frame(links)

  #find the nearest parent
  others<-links$parent_label[!links$parent_label%in%nodes$label]
  others=(others[!duplicated(others)])
  tree%>%filter(label%in%others)->o_nodes
  others=o_nodes$label
  o_nodes=o_nodes$node
  for (i in seq_along(o_nodes)){
    tidytree::ancestor(tree,o_nodes[i])%>%pull(label)->tmp
    links[links$parent_label==others[i],"parent_label"]=rev(tmp)[rev(tmp)%in%nodes$label][1]
  }

  le[le%in%nodes$level]->mytax
  taxRank_to_depth<-setNames(seq_along(mytax)-1, mytax)
  nodes$depth<-taxRank_to_depth[nodes$level%>%as.character()]

  links$IDsource <- match(links$parent_label, nodes$label) - 1
  links$IDtarget <- match(links$label, nodes$label) - 1
  #na.omit(links)->links

  sankeyD3::sankeyNetwork(Links = links, Nodes =nodes ,
                          Source = "IDsource", Target = "IDtarget",Value = "abundant",
                          NodeID = "label",NodeGroup = "label",NodePosX = "depth",NodeValue="abundant",
                          iterations = 1000, xAxisDomain =mytax,align = "none",
                          fontSize = 12,linkGradient = TRUE,
                          nodeWidth = 15,nodeCornerRadius = 5,highlightChildLinks = T,
                          orderByPath = TRUE,scaleNodeBreadthsByString = TRUE,
                          numberFormat = "pavian",dragY = T,nodeShadow = T,
                          doubleclickTogglesChildren = TRUE,width = width,height = height,...)

}

#' Plot a sunburst
#'
#' @return sunburst
#' @export
#' @rdname sangji_plot
#' @examples
#' sunburst(tree)
sunburst<-function(tree){
  tree%>%filter(x<6)->sangji_dat
  #旭日图
  library(plotly)
  fig <- plot_ly(
    #定义所有级别各类的标签
    labels = sangji_dat$label,
    #定义所有级别各类的父级，与上面定义的标签一一对应
    parents = sangji_dat$parent_label,
    #定义各分类的值（一一对应）
    values = sangji_dat$abundant,
    #指定图表类型：sunburst
    type = 'sunburst'
  )
  fig
}



