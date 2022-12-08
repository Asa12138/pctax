#some useful plots

#' Plot correlation
#'
#' @param env dataframe1
#' @param env2 dataframe2 (default:NULL)
#' @param mode plot mode (1~3)
#' @param method one of "pearson","kendall","spearman"
#'
#' @export
#'
#' @examples
#' data(otutab)
#' cor_plot(metadata[,12:19])
#' cor_plot(t(otutab)[,1:50],mode=3)
cor_plot<-function(env,env2=NULL,mode=1,method = "pearson"){
  lib_ps("ggcor")
  library(dplyr)
  set_scale(c("#6D9EC1", "white", "#E46726"),type = "gradient2n")
  if(is.null(env2)){
    if(mode==1){p<-quickcor(env, method = method,cor.test = T) +
      geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
      geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
      geom_abline(slope = -1, intercept=ncol(env)+1)
    detach('package:ggcor')
    return(p)
    }

    if(mode==2) {
      p<-env%>%quickcor(circular = TRUE, cluster = TRUE, open = 45,
                        method = method,cor.test = T) +
        geom_colour(colour = "white", size = 0.125) +
        anno_row_tree() +
        anno_col_tree() +
        set_p_xaxis() +
        set_p_yaxis()
      detach('package:ggcor')
      return(p)
    }

    if(mode==3){
      lib_ps("corrplot")
      ggcor::correlate(env, method = method,cor.test = T,p.adjust = T,p.adjust.method = "fdr")->res2
      rownames(res2$p.value)<-rownames(res2$r);colnames(res2$p.value)<-colnames(res2$r)
      #pls package also has a function called corrplot
      corrplot::corrplot(res2$r, order = "hclust", p.mat = res2$p.value, sig.level = 0.05, insig = "blank",
                         diag = FALSE, tl.cex=0.5, addrect = 17, method="color", outline=TRUE,
                         col=brewer.pal(n=10, name="PuOr"),tl.srt=45, tl.col="black")
    }
  }
  else{
    if(mode==1){p<-quickcor(env,env2, method = method,cor.test = T) +
      geom_square(data = get_data(show.diag = FALSE)) +
      geom_mark(data = get_data(show.diag = FALSE), size = 2.5)
    detach('package:ggcor')
    return(p)
    }

    if(mode==2) {
      p<-quickcor(env,env2,circular = TRUE, cluster = TRUE, open = 45,
                  method = method,cor.test = T) +
        geom_colour(colour = "white", size = 0.125) +
        anno_row_tree() +
        anno_col_tree() +
        set_p_xaxis() +
        set_p_yaxis()
      detach('package:ggcor')
      return(p)
    }
  }
}


#' Plot a stack plot,
#'
#' @param otutab otutab
#' @param metadata metadata
#' @param topN plot how many top species
#' @param groupID one group name of columns of metadata
#' @param shunxu should order the samples by the top1 abundance
#' @param relative transfer to relative or absolute
#' @param style "group" or "sample"
#' @param sorted should legend be sorted by "abunance"
#' @param flow should plot a flow plot?
#' @param others should plot others?
#' @param pmode fill/stack/dodge
#'
#' @export
#'
#' @examples
#' data(otutab)
#' stackplot(otutab,metadata,groupID="Group")
#' stackplot(otutab,metadata,groupID="Id",shunxu=T,flow=T)
#'
#' hclust(dist(t(otutab)))%>%ape::as.phylo()%>%as_tibble()->s_tree
#' full_join(s_tree,metadata,by=c("label"="Id"))->s_tree
#' library(ggtree)
#' ggtree(tidytree::as.treedata(s_tree))+geom_tippoint(aes(col=Group,shape=Group),size=2)->p1
#' stackplot(hebing(otutab,taxonomy$Phylum,1),metadata,groupID ='Id',topN = 10,others = T,flow = T)+
#'   coord_flip()+ scale_y_continuous(expand=c(0,0)) +xlab("")->pp2
#' pp2%>%aplot::insert_left(p1, width=.3)
#'

stackplot<-function (otutab, metadata, topN = 8, groupID = "Group", shunxu=F,relative=T,
                      style = "group", sorted = "abundance",flow=F,others=T,pmode='stack') {
  #用来画物种堆积图，适合处理各种OTU类似数据，输入metatab作为分组依据。style可以选择group或者sample
  #others=T用来选择是否画出除TopN外的其他，pmode可选择fill/stack/dodge
  lib_ps("ggplot2", "reshape2","scales")

  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx, , drop = F]
  otutab = otutab[, rownames(metadata)]

  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  mean_sort = as.data.frame(otutab[(order(-rowSums(otutab))), ])
  if (nrow(mean_sort)>topN){
    other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
    mean_sort = mean_sort[1:(topN - 1), ]
    mean_sort = rbind(mean_sort, other)
    rownames(mean_sort)[topN] = c("Other")
  }

  if (style == "sample") {

    mean_sort$Taxonomy = rownames(mean_sort)
    data_all = as.data.frame(melt(mean_sort, id.vars = c("Taxonomy")))
    if(relative){
      data_all <- data_all  %>%
        group_by(variable, Taxonomy) %>%
        summarise(n = sum(value)) %>%
        mutate(value = n / sum(n))
    }


    if (sorted == "abundance") {
      data_all$Taxonomy = factor(data_all$Taxonomy, levels = rownames(mean_sort))
    }


    data_all = merge(data_all, sampFile, by.x = "variable",
                     by.y = "row.names")
    if (!others){
      data_all<-data_all[data_all$Taxonomy!='Other',]
    }
    if(!flow){
      p = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) +
        geom_bar(stat = "identity",width = 1,position = pmode) +
        facet_grid(~group, as.table = FALSE,
                   switch = "both", scales = "free", space = "free") +
        theme(strip.background = element_blank()) +
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
        xlab(groupID) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    }
    else{
      lib_ps("ggalluvial")
      p = ggplot(data_all, aes(x = variable, y = value,alluvium = Taxonomy, fill = Taxonomy)) +
        ggalluvial::geom_flow(stat="alluvium", lode.guidance = "frontback", color = "darkgray") +
        ggalluvial::geom_stratum(stat="alluvium") +
        facet_grid(~group, as.table = FALSE,
                   switch = "both", scales = "free", space = "free") +
        theme(strip.background = element_blank()) +
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
        xlab(groupID) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    }
    if(relative)p=p+scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance (%)")
    else p=p+ylab("Reads")
    p
  }

  else {
    mat_t = t(mean_sort)
    aggregate(mat_t,by=list(sampFile$group),mean)%>%melt(.,id=1)->data_all
    colnames(data_all)<-c('variable','Taxonomy','value')


    data_all$value = as.numeric(data_all$value)
    as.factor(data_all$variable)->data_all$variable

    if(relative){
      data_all <- data_all  %>%
        group_by(variable, Taxonomy) %>%
        summarise(n = sum(value)) %>%
        mutate(value = n / sum(n))
    }

    if (!others){
      data_all<-data_all[data_all$Taxonomy!='Other',]
    }
    if (sorted == "abundance") {
      data_all$Taxonomy = factor(data_all$Taxonomy, levels = rownames(mean_sort))
    }
    if(shunxu==1)data_all<-mutate(data_all,variable=factor(variable,levels = (data_all%>%filter(Taxonomy==rownames(mean_sort)[1])%>%
                                                                                arrange(value)%>%as.data.frame())[,1]%>%as.character()))
    if(shunxu=='other')data_all<-mutate(data_all,variable=factor(variable,levels = (data_all%>%filter(Taxonomy=='Other')%>%
                                                                                      arrange(value)%>%as.data.frame())[,1]%>%as.character()))
    if(!flow){
      p = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) +
        geom_bar(stat = "identity",position = pmode,width = 0.7) +
        xlab(groupID)
    }
    else{
      lib_ps("ggalluvial")
      p = ggplot(data_all, aes(x = variable, y = value,alluvium = Taxonomy, fill = Taxonomy)) +
        ggalluvial::geom_flow(stat="alluvium", lode.guidance = "frontback", color = "darkgray") +
        ggalluvial::geom_stratum(stat="alluvium") +
        xlab(groupID)
    }

    if(relative)p=p+scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance (%)")
    else p=p+ylab("Reads")
    p
  }

}


#' Heatmap by ggplot
#'
#' @param otutab otutab
#' @param topN plot how many top species
#' @param rowcl cluster the row?
#' @param colcl cluster the column?
#' @param r_ann row annotation
#' @param c_ann column annotation
#'
#' @return a ggplot
#' @export
#'
#' @examples
#' data(otutab)
#' heat_otu(otutab,topN = 30,r_ann = otutab[1:30,1:2],c_ann = metadata[,c(2,4)])
#'
heat_otu<-function(otutab,topN=20,rowcl=T,colcl=T,r_ann=NULL,c_ann=NULL){
  lib_ps("ggplot2","gplots","ggnewscale","aplot","tidyr","ggtree","ape")
  otutab[1:topN,]->d
  rownames(d)->d$otu
  dd<-tidyr::gather(d,1:ncol(d)-1, key="sample", value='value')

  p <- ggplot(dd, aes(sample,otu, fill=value)) + geom_tile() +
    scale_fill_gradientn(colours = rev(gplots::redblue(50))) +
    scale_y_discrete(position="right") +
    theme_minimal() +
    xlab(NULL) + ylab(NULL)

  if(!is.null(r_ann)){
    ca1 <- r_ann
    rownames(ca1)->ca1$Id
    pc1=ggplot() + geom_tile(data = ca1, aes(x=colnames(ca1)[1], y=Id, fill=get(colnames(ca1)[1])))+
      labs(fill = colnames(ca1)[1])
    if(ncol(ca1)>2){
      for (i in 2:(ncol(ca1)-1)){
        pc1=pc1+ggnewscale::new_scale_fill()+
          geom_tile(data = ca1, aes(y=Id,x=colnames(ca1)[i], fill=get(colnames(ca1)[i])))+
          labs(fill = colnames(ca1)[i])
      }
    }
    pc1=pc1+scale_y_discrete(position="right") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      xlab(NULL) + ylab(NULL)
    p=p %>% insert_left(pc1, width=.1)

  }
  if(rowcl){
    hclust(dist(otutab[1:topN,]))%>%ape::as.phylo()->a
    p=p %>% insert_left(ggtree(a,branch.length = "none"), width=.1)
  }
  if(!is.null(c_ann)){
    ca <- c_ann
    rownames(ca)->ca$Id

    pc=ggplot() + geom_tile(data = ca, aes(x=Id,y=colnames(ca)[1], fill=get(colnames(ca)[1])))+
      labs(fill = colnames(ca)[1])

    if(ncol(ca)>2){
      for (i in 2:(ncol(ca)-1)){
        pc=pc+ggnewscale::new_scale_fill()+
          geom_tile(data = ca, aes(x=Id,y=colnames(ca)[i], fill=get(colnames(ca)[i])))+
          labs(fill = colnames(ca)[i])
      }
    }
    pc=pc+scale_y_discrete(position="right") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      xlab(NULL) + ylab(NULL)

    p=p %>% insert_top(pc, height=.1)
  }
  if(colcl){
    hclust(dist(t(otutab[1:topN,])))%>%ape::as.phylo()->b
    p=p %>% insert_top(ggtree(b,branch.length = "none") + layout_dendrogram(), height=.1)
    }
  return(p)
}


#' Areaplot
#'
#' @export
#' @rdname stackplot
#' @examples
#'areaplot(otutab,metadata,topN = 10,others = F)
areaplot<-function (otutab, metadata, topN = 8, groupID = "Date",relative=T,
                    sorted = "abundance",others=T) {
  #用来画物种堆积图，适合展示时间数据，输入metatab作为分组依据。style可以选择group或者sample
  #others=T用来选择是否画出除TopN外的其他，pmode可选择fill/stack/dodge
  lib_ps("ggplot2", "reshape2","scales")

  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx, , drop = F]
  otutab = otutab[, rownames(metadata)]

  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  mean_sort = as.data.frame(otutab[(order(-rowSums(otutab))), ])
  if (nrow(mean_sort)>topN){
    other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
    mean_sort = mean_sort[1:(topN - 1), ]
    mean_sort = rbind(mean_sort, other)
    rownames(mean_sort)[topN] = c("Other")
  }

  {
    mat_t = t(mean_sort)
    aggregate(mat_t,by=list(sampFile$group),mean)%>%melt(.,id=1)->data_all
    colnames(data_all)<-c('variable','Taxonomy','value')

    if (sorted == "abundance") {
      data_all$Taxonomy = factor(data_all$Taxonomy, levels = rownames(mean_sort))
    }

    data_all$value = as.numeric(data_all$value)
    data_all$variable<-as.Date(data_all$variable)

    if(relative){
      data_all <- data_all  %>%
        group_by(variable, Taxonomy) %>%
        summarise(n = sum(value)) %>%
        mutate(value = n / sum(n))
    }

    datebreaks <- unique(data_all$variable)

    if (!others){
      data_all<-data_all[data_all$Taxonomy!='Other',]
    }

    p=ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) +
      geom_area() +scale_x_date(breaks = datebreaks)+ ylab(label = 'Relative Abundance(%)')+xlab(NULL)+
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
            plot.margin = margin(t = 0.5,  r = 1,  b = 0.5, l = 0.5,unit = "cm"))
    p
  }
}

#' Circle plot
#'
#' @param rowcol colors codes
#' @param colcol colors codes
#' @rdname stackplot
#' @return plot
#' @export
#'
#' @examples
#' tax_circlize(otutab,metadata,groupID ='Site',topN = 10,others = F)
tax_circlize<-function (otutab, metadata,topN = 5, groupID = "Group",others=T,rowcol=c(),colcol=c()) {
  lib_ps("ggplot2", "reshape2", "circlize","RColorBrewer")
  #取出metatab和otutab有重叠的部分
  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx, , drop = F]
  otutab = otutab[, rownames(metadata)]
  #生成一列的dataframe，只有一个分组信息
  sampFile = as.data.frame(metadata[, groupID]%>%as.factor(), row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  mean_sort = as.data.frame(otutab[(order(-rowSums(otutab))), ])
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  #可能有的画图数据全部要展示，不要others，也不是相对丰度
  mean_sort2 = mean_sort[1:topN, ]

  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort, other)
  rownames(mean_sort)[topN] = c("Other")

  if (!others){
    mean_sort = mean_sort2}

  mat_t = t(mean_sort)
  mat_t2 = merge(sampFile, mat_t, by = "row.names")
  mat_t2 = mat_t2[, c(-1)]
  mat_mean = aggregate(mat_t2[, -1], by = mat_t2[1], FUN = mean)
  df = t(mat_mean[, -1])
  geno = mat_mean$group
  colnames(df) = mat_mean$group

  #设置颜色这样做
  grid.col = NULL

  if (length(rowcol)==0) grid.col[rownames(df)] = brewer.pal(dim(df)[1], "Set3")
  else grid.col[rownames(df)] =rowcol

  if (length(colcol)==0) grid.col[colnames(df)] = brewer.pal(dim(df)[2], "Accent")
  else grid.col[colnames(df)] =colcol

  chordDiagram(df, directional = TRUE, diffHeight = 0.03, grid.col = grid.col, transparency = 0.5)
  legend("left", pch = 20, legend = rownames(df), col = grid.col[rownames(df)],
         bty = "n", cex = 1, pt.cex = 3, border = "black")
  legend("right", pch = 20, legend = colnames(df), col = grid.col[colnames(df)],
         bty = "n", cex = 1, pt.cex = 3, border = "black")
  detach('package:circlize')
}

#' Pie plot
#'
#' @param otutab otutab
#' @param n topn
#'
#' @return ggplot
#' @export
#'
#' @examples
#'tax_pie(otutab,n = 7)
tax_pie<-function(otutab,n=6){
  lib_ps("RColorBrewer","ggpbur")
  rowSums(otutab)->a
  if(length(a)>n){
    sort(a,decreasing = T)[1:n-1]->b
    other=sum(sort(a,decreasing = T)[n:length(a)])
    b<-c(b,other)
    names(b)[length(b)]<-'Others'}
  else b<-a
  myPalette <- brewer.pal(n, "Paired")
  # You can change the border of each area with the classical parameters:
  # pie(b , labels = paste0(names(b),"\n(",round(b/sum(b)*100,2),"%)"), border="white",
  #     col=myPalette,radius = 1,main = main)
  df=data.frame(va=b,labels = paste0(names(b),"\n(",round(b/sum(b)*100,2),"%)"))
  ggpie(df,'va',fill=myPalette,label = "labels",grepl=T)
}

#' Radar plot
#'
#' @param otu_time
#'
#' @export
#'
#' @examples
#' tax_radar(otu_time)
tax_radar<-function(otu_time){
  lib_ps("ggradar","scales")
  otu_time[1:4,]%>%
    mutate_all(scales::rescale) %>%cbind(tax=rownames(.),.)%>%
    ggradar::ggradar()
}

#' Word cloud
#'
#' @param aa
#'
#' @export
#'
#' @examples
#'tax_wordcloud(taxonomy$Genus)
tax_wordcloud<-function(aa){
  lib_ps("pcutils","wordcloud2")
  remove_unclassfied<-\ (taxdf) {
    taxdf[grepl.data.frame("Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig",
                           taxdf, ignore.case = TRUE)] <- NA
    return(taxdf)
  }
  sort(table(aa),decreasing = TRUE)[1:50]%>%as.data.frame()%>%
    remove_unclassfied()%>%na.omit()%>%wordcloud2::wordcloud2(.,size=.7)
}

#' Triangle plot
#'
#' @param otutab otutab
#' @param group group
#' @param scale default:F
#' @param class
#'
#' @export
#'
#' @examples
#' data(otutab)
#'triangp(otutab,metadata$Group,class=taxonomy$Phylum)+theme_classic()
triangp<-function(otutab,group,scale=F,class=NULL){
  lib_ps("ggtern","vegan")
  group%>%as.factor()->group
  if (nlevels(group)!=3)stop("group is not 3, can't plot trip")
  hebing(otutab,group,act = 'mean')->tmp
  if (scale){vegan::decostand(tmp,'hellinger',2)->tmp}
  tmp%>%as.data.frame()%>%mutate(sum=rowSums(.))->tmp1
  colnames(tmp1)[1:3]<-c('KO','OE','WT')
  if (is.null(class)){
    p=ggtern(tmp1,aes(x=KO,y=OE,z=WT)) +
      geom_point(aes(size=sum,col=class))+#define data geometry
      theme_showarrows() +labs(x=names(tmp)[1],y=names(tmp)[2],z=names(tmp)[3])
    return(p)
  }
  else {
    tmp1$class =class
    p=ggtern(tmp1,aes(x=KO,y=OE,z=WT)) +
      geom_point(aes(size=sum,col=class))+#define data geometry
      theme_showarrows() +labs(x=names(tmp)[1],y=names(tmp)[2],z=names(tmp)[3])
    return(p)
  }
}

