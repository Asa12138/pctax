#phylogenic tree==============
#' Complete a taxonomy table
#'
#' @param taxdf taxonomy table
#' @param type species
#' @import dplyr tibble
#' @return a good taxonomy table
#' @export
#'
#' @examples
#'taxmat = matrix(sample("onelevel", 7*2, replace = TRUE), nrow = 2, ncol = 7)%>%as.data.frame()
#'colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#'fillNAtax(taxmat)

fillNAtax<-function (taxdf, type = "species") {
  lib_ps("dplyr","tibble")
  taxaclass=c( "Kingdom","Phylum","Class","Order","Family","Genus","Speies" ,"Rank8", "Rank9" ,"Rank10")

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
  repduplicatedtaxcheck<-function (taxdf) {
    for (i in seq_len(7)) {
      taxdf <- duplicatedtaxcheck(taxdf) %>% column_to_rownames(var = "rowname")
    }
    return(taxdf)
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

#' From a dataframe to construct a phylo
#'
#' @param taxa dataframe
#'
#' @return phylo
#' @export
#'
#' @examples
#' data(otutab)
#' makeNewick(taxonomy)->spe_tree
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
#' @param otutab otutab, rowname==colname(taxonomy)
#' @param level 1~7
#'
#' @return a treedata
#' @export
#'
#' @examples
#' data(otutab)
#' ann_tree(taxonomy,otutab)->tree
#' easy_tree(tree)
ann_tree<-function(f_tax,otutab,level=7){
  if(any(rownames(f_tax)!=rownames(otutab)))stop("rowname not match")
  otutab%>%rowSums()->num
  res=data.frame(label="",abundant=sum(num))
  for (i in 1:level){
    aggregate(num,by=list(f_tax[,i]),sum)->tmp
    colnames(tmp)<-colnames(res)
    res=rbind(res,tmp)
  }
  le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  makeNewick(f_tax[,1:level])%>%fortify()->tree
  tree$level<-le[ceiling(tree$branch)+1]
  #tree$level<-factor(tree$level,levels = le)
  left_join(tree,tree[,c("node","label")],by=c("parent"="node"))%>%dplyr::rename(label=label.x,parent_l=label.y)->tree1
  left_join(tree1,res)->tree2
  left_join(tree2,f_tax[,1:level]%>%distinct(),by=c("label"=colnames(f_tax)[level]))->tree3
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
easy_tree<-function(tree){
  if(!requireNamespace("treedataverse"))BiocManager::install("YuLab-SMU/treedataverse");
  library(treedataverse)
  #可以剪切部分树
  tree%>%mutate(label1=sub(".?__","",label))->tree2
  ggtree(tree2,layout = 'radial',size=0.5)+
    geom_tiplab(color="black",size=1.5,offset = 1, show.legend = FALSE)+
    geom_highlight(data = subset(tree2,branch==1.5),aes(node = node, fill = label1),alpha=0.5)+
    ggnewscale::new_scale_fill()+
    #scale_fill_d3()+
    geom_fruit(
      geom=geom_tile,
      mapping=aes(y=label, fill=log(abundant)),
      stat="identity",offset = 0.1
    )+scale_fill_gradient2(low = "blue",mid = "#FFFFFF",high = "red")
}

#' Plot a sankey
#'
#' @param tree ann_tree result
#' @param top_N each level has top_N
#'
#' @export
#'
#' @examples
#' sangji_plot(tree)
sangji_plot<-function(tree,top_N=5){
  if(!requireNamespace("sankeyD3"))devtools::install_github("fbreitwieser/sankeyD3");
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
  tree%>%group_by(level)%>%top_n(top_N,abundant)%>%ungroup()->sangji_dat
  sangji_dat%>%filter(label!="")%>%select(label,level,abundant)->nodes
  le[le%in%nodes$level]->mytax
  taxRank_to_depth<-setNames(seq_along(mytax)-1, mytax)
  nodes$depth<-taxRank_to_depth[nodes$level%>%as.character()]

  sangji_dat%>%filter(parent_l!="")%>%select(label,parent_l,abundant)->links
  others<-links$parent_l[!links$parent_l%in%nodes$label]
  tree%>%filter(label%in%others)%>%pull(node)->o_nodes
  for (i in seq_along(o_nodes)){
    ancestor(tree,o_nodes[i])%>%pull(label)->tmp
    links[links$parent_l==others[i],"parent_l"]=rev(tmp)[rev(tmp)%in%nodes$label][1]
  }

  links$IDsource <- match(links$parent_l, nodes$label) - 1
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
                          doubleclickTogglesChildren = TRUE,width = 2000,height = 1000)

}

#' Plot a sunburst
#'
#' @return
#' @export
#' @rdname sangji_plot
#' @examples
#' sunburst(tree,10)
sunburst<-function(tree,top_N=5){
  le=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tree%>%group_by(level)%>%top_n(top_N,abundant)%>%ungroup()->sangji_dat
  sangji_dat%>%filter(label!="")%>%select(label,level,abundant)->nodes
  le[le%in%nodes$level]->mytax
  taxRank_to_depth<-setNames(seq_along(mytax)-1, mytax)
  nodes$depth<-taxRank_to_depth[nodes$level%>%as.character()]

  sangji_dat%>%filter(parent_l!="")%>%select(label,parent_l,abundant)->links
  links$IDsource <- match(links$parent_l, nodes$label) - 1
  links$IDtarget <- match(links$label, nodes$label) - 1
  na.omit(links)->links
  #旭日图
  library(plotly)
  fig <- plot_ly(
    #定义所有级别各类的标签
    labels = links$label,
    #定义所有级别各类的父级，与上面定义的标签一一对应
    parents = links$parent_l,
    #定义各分类的值（一一对应）
    values = links$abundant,
    #指定图表类型：sunburst
    type = 'sunburst'
  )
  fig
}
