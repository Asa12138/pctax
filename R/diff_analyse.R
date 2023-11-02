#different analyse==============

#' Difference analysis
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param ctrl the control group, one level of groups
#' @param method one of "deseq2","edger","limma","t.test","wilcox.test"
#' @param log do log transfer for limma?
#' @param add_mini add_mini when calculate the logFC. e.g (10+0.1)/(0+0.1), default `0.5*min(abundance)`
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' data(otutab,package = "pcutils")
#' diff_da(otutab,metadata["Group"],method="deseq2")->res
#' volcano_p(res)
#' volcano_p(res,mode=2)
diff_da<-function(otutab,group_df,ctrl=NULL,method="deseq2",log=T,add_mini=NULL){
  if(length(method)>1){
    all_res=list()
    for (i in method) {
      pcutils::dabiao("Try ",i)
      res=tryCatch({
        diff_da(otutab,group_df,ctrl,method = i,log,add_mini)},
        error = function(e) {
        NULL
      })
      all_res[[i]]=res
    }
    return(all_res)
  }
  method=match.arg(method,c("deseq2","edger","limma","t.test","wilcox.test"))
  idx = rownames(group_df) %in% colnames(otutab)
  group_df = group_df[idx, , drop = F]
  otutab = otutab[, rownames(group_df),drop=F]
  group_df%>%dplyr::rename(Group=1)->meta
  meta$Group=factor(meta$Group)
  if(!is.null(ctrl))meta$Group<-relevel(meta$Group,ctrl)
  else ctrl=levels(meta$Group)[1]

  if(method=="deseq2"){
    lib_ps("DESeq2",library = F)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = otutab, colData = meta, design = ~Group)#构建 DESeqDataSet 对象
    dds <- DESeq2::DESeq(dds) #差异分析

    #results(dds,name = resultsNames(dds)[2])%>%as.data.frame()
    res=data.frame()
    for (i in 2:length(DESeq2::resultsNames(dds))) {
      DESeq2::results(dds,name = DESeq2::resultsNames(dds)[i])%>%as.data.frame()->tmp
      tibble::rownames_to_column(tmp,"tax")->tmp
      res=rbind(res,data.frame(tmp,compare=DESeq2::resultsNames(dds)[i]))
    }
    res$compare=gsub("Group_","",res$compare)%>%gsub("_vs_","-",.)
    res$method="DESeq2"
  }

  if(method=="edger"){
    lib_ps("edgeR",library = F)
    res=data.frame()

    for (i in 2:nlevels(meta$Group)){
      rbind(filter(meta,Group==levels(meta$Group)[i]),
            filter(meta,Group==levels(meta$Group)[1]))->meta1
      group=meta1$Group
      otutab_f=otutab[,rownames(meta1)]

      #数据预处理
      #（1）构建 DGEList 对象
      dgelist <- edgeR::DGEList(counts = otutab_f, group = group)
      #（2）过滤 low count 数据，例如 CPM 标准化（推荐）
      # keep <- rowSums(cpm(dgelist) > 1 ) >= 2
      # dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
      #（3）标准化，以 TMM 标准化为例
      dgelist_norm <- edgeR::calcNormFactors(dgelist, method = 'TMM')
      design <- stats::model.matrix(~factor(meta1$Group))

      #（1）估算基因表达值的离散度
      dge <- edgeR::estimateDisp(dgelist_norm, design, robust = TRUE)
      #（2）模型拟合，edgeR 提供了多种拟合算法
      #负二项广义对数线性模型
      fit <- edgeR::glmFit(dge, design, robust = TRUE)
      lrt <- edgeR::topTags(edgeR::glmLRT(fit), n = nrow(dgelist$counts))
      lrt$table->tmp
      tibble::rownames_to_column(tmp,"tax")->tmp
      res=rbind(res,data.frame(tmp,compare=paste(levels(meta$Group)[c(i,1)],collapse = "-")))
    }
    colnames(res)=c("tax","log2FoldChange","logCPM","LR","pvalue","padj","compare")
    res$method="edgeR"
  }

  if(method=="limma"){
    #limma
    lib_ps("limma",library = F)
    otutab_log=otutab
    if(log)otutab_log=log(otutab+1)
    #log后的数据哦

    #df.matrix <- makeContrasts(KO - WT , levels = list)
    res=data.frame()
    for (i in 2:nlevels(meta$Group)){
      rbind(filter(meta,Group==levels(meta$Group)[i]),
            filter(meta,Group==levels(meta$Group)[1]))->meta1
      group=meta1$Group
      otutab_f=otutab_log[,rownames(meta1)]

      list <- stats::model.matrix(~0+factor(meta1$Group))  #把group设置成一个model matrix
      colnames(list) <- c("c","t")
      df.fit <- limma::lmFit(otutab_f, list)  ## 数据与list进行匹配

      df.matrix <-  limma::makeContrasts(t-c,levels = list)

      fit <-  limma::contrasts.fit(df.fit, df.matrix)
      fit <-  limma::eBayes(fit)
      tempOutput <-  limma::topTable(fit,n = Inf,coef = 1)
      tibble::rownames_to_column(tempOutput,"tax")->tempOutput
      res=rbind(res,data.frame(tempOutput,compare=paste(levels(meta$Group)[c(i,1)],collapse = "-")))
    }
    colnames(res)=c("tax","log2FoldChange","AveExpr","t","pvalue","padj","B","compare")
    res$method="limma"
  }

  if(method%in%c("t.test","wilcox.test")){
    lib_ps("ggpubr",library = F)
    t(otutab)%>%data.frame(.,check.names = F)%>%dplyr::mutate(Group=meta$Group)->dat
    reshape2::melt(dat,id.vars = "Group")->dat
    ggpubr::compare_means(formula = value ~ Group, data = dat,
                          group.by = "variable", method = method,
                          p.adjust.method = "fdr",ref.group = ctrl)->x.all
    x.all=dplyr::arrange(x.all,group1,group2)
    x.all%>%rename(tax="variable","padj"="p.adj")->x.all
    x.all$tax=as.character(x.all$tax)
    x.all$p.format=as.numeric(x.all$p.format)
    x.all=dplyr::filter(x.all,group1==ctrl)
    x.all$compare=paste0(x.all$group2,"-",x.all$group1)

    pcutils::hebing(otutab,meta$Group)->tmp
    if(is.null(add_mini))add_mini=min(tmp[tmp>0])*0.5
    logfc=data.frame()
    for (i in 2:nlevels(meta$Group)){
      logfc=rbind(logfc,data.frame(tax=rownames(tmp),
                 compare=paste0(levels(meta$Group)[i],"-",levels(meta$Group)[1]),
                 log2FoldChange=log2((tmp[,levels(meta$Group)[i]]+add_mini)/(tmp[,levels(meta$Group)[1]]+add_mini))
      ))
    }
    res=dplyr::left_join(x.all,logfc)
  }

  res[which(res$log2FoldChange >=1  & res$padj < 0.05),'sig'] <- 'up'
  res[which(res$log2FoldChange <= -1 & res$padj < 0.05),'sig'] <- 'down'
  res[which(is.na(res$sig)),'sig'] <- 'none'
  res%>%dplyr::mutate(tax1=ifelse(sig%in%c("up","down"),tax,""))->res

  return(res=res)
}

#' Volcano plot for difference analysis
#'
#' @param res result of `diff_da` which have colnames: tax, log2FoldChange, padj, compare, sig
#' @param logfc log_fold_change threshold
#' @param adjp adjust_p_value threshold
#' @param mode 1:normal; 2:multi_contrast
#' @param number show the tax number
#' @param text text, T
#'
#' @return ggplot
#' @export
#'
#' @seealso \code{\link{diff_da}}
volcano_p<-function(res,logfc=1,adjp=0.05,text=T,mode=1,number=F){
  lib_ps("ggrepel","reshape2",library = F)

  this_method=unique(res$method)
  if(length(this_method)!=1)stop("Wrong method column")

  res$sig=NULL
  res[which(res$log2FoldChange >=logfc  & res$padj < adjp),'sig'] <- 'up'
  res[which(res$log2FoldChange <= -logfc & res$padj < adjp),'sig'] <- 'down'
  res[which(is.na(res$sig)),'sig'] <- 'none'
  res%>%dplyr::mutate(tax1=ifelse(sig%in%c("up","down"),tax,""))->res
  res=dplyr::arrange(res,-log2FoldChange)
  res$label <- ifelse(res$sig%in%c("up","down"),"Sig","Non-sig")
  cols=setNames(c("#2f5688","#BBBBBB","#CC0000"),c("down","none","up"))

  if(number){
    label=dplyr::count(res,sig)%>%dplyr::mutate(label=paste0(sig,": ",n))
    res$sig=setNames(label$label,label$sig)[res$sig]
    cols=setNames(cols,setNames(label$label,label$sig)[names(cols)])
  }

  if(mode==1){
    #unique(res$compare)
    #一对比较的火山图
    res%>%dplyr::filter(!is.na(padj))->dat

    pp<-ggplot(dat,aes(x=log2FoldChange,y=-log10(padj),color=sig))+
      geom_point()+
      scale_color_manual(values=cols,na.value ="#BBBBBB")+  #确定点的颜色
      facet_wrap(.~compare,scales = "free")+labs(title = this_method)+
      theme_bw()+  #修改图片背景
      theme(
        legend.title = element_blank()  #不显示图例标题
      )+
      ylab('-log10 (p-adj)')+  #修改y轴名称
      xlab('log2 (FoldChange)')+  #修改x轴名称
      geom_vline(xintercept=c(-logfc,logfc),lty=3,col="black",lwd=0.5) +  #添加垂直阈值|FoldChange|>2
      geom_hline(yintercept = -log10(adjp),lty=3,col="black",lwd=0.5)  #添加水平阈值padj<0.05

    if(text)pp=pp+
      ggrepel::geom_text_repel(aes(label=tax1),size=2)
  }
  if(mode==2){
    #多对比较的火山图
    res%>%dplyr::filter(abs(log2FoldChange)>logfc)%>%dplyr::filter(!is.na(padj))->dat
    res%>%dplyr::group_by(compare)%>%dplyr::summarise(max(log2FoldChange))%>%reshape2::melt()->bardf
    res%>%dplyr::group_by(compare)%>%dplyr::summarise(min(log2FoldChange))%>%reshape2::melt()->bardf1

    p1 <- ggplot()+
      geom_bar(data = bardf,
               mapping = aes(x = compare,y = value),stat ="identity",
               fill = "#dcdcdc",alpha = 0.6)+
      geom_bar(data = bardf1,
               mapping = aes(x = compare,y = value),stat ="identity",
               fill = "#dcdcdc",alpha = 0.6)
    p2<-p1+geom_jitter(data = dat,
                       aes(x = compare, y = log2FoldChange, color = label),
                       size = 1.2,
                       width =0.4)

    if(text)p2=p2+
      ggrepel::geom_text_repel(data=filter(dat,tax1!=""),aes(x = compare, y = log2FoldChange, label=tax1),col="red",size=3,
                               force = 1.2,arrow = arrow(length = unit(0.008, "npc"),
                                                         type = "open", ends = "last"))

    coldf<-data.frame(x=unique(res$compare),y=0)
    p3<-p2+geom_tile(data = coldf,
                     aes(x=x,y=y,fill=x),
                     height=0.4,
                     color = "black",
                     alpha = 0.6,
                     show.legend = F)+
      labs(x="Compares",y="log2 (FoldChange)")+
      geom_text(data=coldf,
                aes(x=x,y=y,label=x),
                size =6,
                color ="white")+
      scale_color_manual(name=NULL,
                         values = c("black","red"))


    pp=p3+theme_minimal()+
      #coord_cartesian(ylim = c(-2,4))+
      theme(
        axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        #legend.direction = "h",
        legend.text = element_text(size = 15)
      )
  }
  return(pp)
}

#' Difference analysis
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' data(otutab,package = "pcutils")
#' multi_bar(otutab[1:10,],metadata["Group"])
multi_bar=function(otutab,group_df,mode=1,text_df=NULL,text_x=NULL,text_angle=-90){
  idx = rownames(group_df) %in% colnames(otutab)
  group_df = group_df[idx, , drop = F]
  otutab = otutab[, rownames(group_df),drop=F]
  group_df%>%dplyr::rename(Group=1)->meta
  #meta$Group=factor(meta$Group)

  t(otutab)%>%data.frame(.,check.names = F)%>%mutate(Group=meta$Group)->dat
  reshape2::melt(dat,id.vars = "Group",variable.name = "tax",value.name = "abundance")->dat
  dat$tax=factor(dat$tax,levels = rownames(otutab))

  low=min(dat$abundance);high=max(dat$abundance)

  if(mode==1){
    meandf=group_by(dat,Group,tax)%>%summarise(mean=mean(abundance),sd=sd(abundance),se=sd(abundance)/sqrt(n()-1))
    p=ggplot(data = meandf,aes(y=tax,x=mean))+
      geom_errorbar(aes(xmin=mean-se,xmax=mean+se,Group=Group),stat = "identity",
                    position = position_dodge(width=0.9),width=0.25)+
      geom_bar(aes(fill=Group),stat = "identity",position = position_dodge(width=0.9),
               width = 0.8,size=0.2,color="black")
    if(is.null(text_x))text_x=-0.05*(high)
  }
  if(mode==2){
    p=ggplot(data = dat,aes(y=tax,x=abundance))+
      geom_boxplot(aes(color=Group),outlier.shape = NA)
    if(is.null(text_x))text_x=low-0.05*(high-low)
  }
  if(!is.null(text_df)){
    text_df=text_df[rownames(otutab),,drop=F]
    text_df$tax=rownames(text_df)
    colnames(text_df)[1]="label"
    p=p+geom_text(data = text_df,mapping = aes(x=text_x,y=tax,label=label),angle=text_angle)
  }
  p
  # lib_ps("ggpubr")
  # lib_ps("ggfun")
  # ggpubr::ggbarplot(dat,y="abundance",x="tax",fill = "Group",add = "mean_se",position = position_dodge(width = 0.8))+coord_flip()+
  #   ggfun::theme_stamp(axis = "x")
}


#' Get mean and type
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @export
get_diff_type=function(otutab,group_df){
  idx = rownames(group_df) %in% colnames(otutab)
  group_df = group_df[idx, , drop = F]
  otutab = otutab[, rownames(group_df),drop=F]

  group=group_df[[1]]%>%as.factor()
  pcutils::hebing(otutab,group)->tmp
  apply(tmp, 1, function(a)which(a==max(a))[[1]])->tmp$group
  apply(tmp, 1, function(a){return(colnames(tmp)[a[[ncol(tmp)]]]) })->tmp$Type
  tmp$Type<-factor(tmp$Type,levels = levels(group))
  tmp$tax=rownames(tmp)
  tmp
}


#' KW test
#'
#' @param otutab otutab
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param method "kruskal.test", see \code{\link[ggpubr]{compare_means}}
#'
#' @return res
#' @export
#'
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' kwtest(otutab,metadata["Group"])->res
#' bbtt(res,pvalue = "p.format")
#' }}
kwtest<-function(otutab,group_df,method = "kruskal.test"){
  lib_ps("ggpubr","reshape2",library = F)
  group=group_df[[1]]%>%as.factor()
  t(otutab)%>%data.frame(.,check.names = F)%>%dplyr::mutate(Group=group)->dat
  reshape2::melt(dat,id.vars = "Group")->dat
  ggpubr::compare_means(formula = value ~ Group, data = dat,  group.by = "variable",
                        method = method, p.adjust.method = "fdr")->x.all
  x.all%>%rename(tax="variable")->x.all
  x.all$tax=as.character(x.all$tax)
  x.all$p.format=as.numeric(x.all$p.format)

  tmp=get_diff_type(otutab,group_df)
  x.all=dplyr::left_join(x.all,tmp,by="tax")
  return(x.all)
}

#' ALDEX
#'
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param otutab otutab
#'
#' @return diff
#' @export
#' @references https://cloud.tencent.com/developer/article/1621879
#'
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' ALDEX(otutab,metadata["Group"])->res
#' res%>%dplyr::top_n(9,-glm.eBH)%>%.[,"tax"]->sig
#' data.frame(t(otutab[sig,]))%>%pcutils::group_box(.,"Group",metadata)
#' }}
ALDEX<-function(otutab,group_df){
  lib_ps("ALDEx2",library = F)
  group=group_df[[1]]%>%as.factor()
  if(nlevels(factor(group))==2){
    x.all <- ALDEx2::aldex(otutab, group, mc.samples=16, test="t", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
    data.frame(tax=rownames(x.all),x.all)->x.all
    x.all%>%mutate(p.signif=ifelse(we.eBH>=0.05,'ns',ifelse(we.eBH>=0.01,'*',ifelse(we.eBH>0.001,'**','***'))))->x.all
  }
  if(nlevels(factor(group))>2){
    x.all <- ALDEx2::aldex(otutab, group, mc.samples=16, test="kw", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
    data.frame(tax=rownames(x.all),x.all)->x.all
    x.all%>%mutate(p.signif=ifelse(glm.eBH>=0.05,'ns',ifelse(glm.eBH>=0.01,'*',ifelse(glm.eBH>0.001,'**','***'))))->x.all
  }

  tmp=get_diff_type(otutab,group_df)
  x.all=dplyr::left_join(x.all,tmp,by="tax")
  return(x.all)
}

# ANCOM2<-function(otutab,group_df){
#   res = ANCOM(otutab, group_df, struc_zero = NULL,main_var = colnames(group_df))
#   res$out->out
#   rownames(out)<-out$taxa_id
#   hebing(otutab,metadata[,group])->tmp
#   apply(tmp, 1, function(a)which(a==max(a))[[1]])->tmp$group
#   apply(tmp, 1, function(a){return(colnames(tmp)[a[[4]]]) })->tmp$Type
#   rename(out,tax=taxa_id)->out
#   cbind(out,tmp[out$tax,])->out
#   cbind(out,otutab[out$tax,])->out
#
#   # Step 3: Volcano Plot
#   n_taxa = nrow(otutab)
#   # Cutoff values for declaring differentially abundant taxa
#   cut_off = c(0.9 * (n_taxa - 1), 0.8 * (n_taxa - 1), 0.7 * (n_taxa - 1), 0.6 * (n_taxa - 1))
#   names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")
#   # Annotation data
#   dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
#   res$fig = res$fig +
#     geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
#     geom_text(data = dat_ann, aes(x = x, y = y, label = label),
#               size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
#   return(res)
# }


#' ggdotchart for diff analysis
#'
#' @param res result of ALDEX or kwtest
#' @param pvalue the name of pvaule
#' @param topN topN
#'
#' @return ggplot
#' @export
bbtt<-function(res,pvalue="glm.eBH",topN=20){
  res%>%dplyr::arrange(get(pvalue))%>%head(topN)%>%
    dplyr::mutate(group=Type,tax=factor(tax,levels = rev(tax)))->pres1
  p1<-ggpubr::ggdotchart(pres1, x = "tax", y = pvalue,color = "group",
                 dot.size = 7,  #palette = classcol, # 修改颜色
                 sorting = "ascending",
                 add = "segments",
                 add.params = list(color = "lightgray", size = 1.5),# 添加棒子
                 label = "Type",# 添加label
                 font.label = list(color = "white", size = 8,vjust = 0.5),
                 rotate = T,
                 ggtheme = pctax_theme# 改变主题
  )+labs(x='')
  p1
  # pres1[,levels(pres1$Type)]%>%pcutils::trans(.,1,method = 'normalize')%>%
  #   dplyr::mutate(tax=rownames(.),tax=factor(tax,levels = rev(tax)))->tmp
  # p2<-reshape2::melt(tmp)%>%
  #   ggplot(.,aes(x=variable,y=tax,fill=value))+geom_raster()+
  #   scale_fill_gradientn(colours = c("#1789E6", "#FFFFFF", "#B3192B"))+
  #   theme_void()+
  #   theme(axis.text.x = element_text(),legend.position = 'top')
  # lib_ps("patchwork",library = F)
  # p2+p1+plot_layout(widths = c(1.5,3))
}

#' RandomForest
#'
#' @param otutab otutab
#' @param topN default: 10
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return diff
#' @export
#'
#' @examples
#' data(otutab,package = "pcutils")
#' suijisenlin(otutab,metadata["Group"])->rf_res
suijisenlin<-function(otutab,group_df,topN=10){
  group=group_df[[1]]

  lib_ps("randomForest",library = F)
  t(otutab)%>%data.frame()%>%mutate(Group=group)%>%
    randomForest::randomForest(Group ~.,data = .,importance=TRUE ,proximity=TRUE)->randomforest

  #plot(randomforest)
  #查看变量的重要性
  randomforest$importance%>%data.frame()%>%tibble::rownames_to_column("tax")->im

  im%>%dplyr::top_n(topN,wt =MeanDecreaseAccuracy)%>%dplyr::arrange(-MeanDecreaseAccuracy)%>%
    dplyr::mutate(tax=factor(tax,levels = rev(tax)))->mim

  imp<-ggpubr::ggscatter(mim,x='MeanDecreaseAccuracy',y='tax',color = '#7175E3',
                size="MeanDecreaseAccuracy",ggtheme = pctax_theme)+ylab('')+
    theme(legend.position = "none")

  # pcutils::hebing(otutab,group =group)->tmp
  # tmp[as.character(mim$tax),]%>%decostand(.,1,method = 'normalize')%>%
  #   mutate(tax=rownames(.),tax=factor(tax,levels = rev(tax)))->tmp
  #
  # p2<-melt(tmp)%>%
  #   ggplot(.,aes(x=variable,y=tax,fill=value))+geom_raster()+
  #   scale_fill_gradientn(colours = c("#003366","white","#990033"))+theme_pubr(legend = 'right')+
  #   theme_void()+
  #   theme(axis.text.x = element_text(),legend.position = 'right')
  # imp<-imp+p2+plot_layout(widths = c(2.5,1.5))
  return(list(im=im,imp=imp))
}

#' Time series analysis
#'
#' @param otu_time otutab hebing by a time variable
#' @param n_cluster number of clusters
#' @param min.std min.std
#'
#' @return time_cm
#' @export
#'
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' otu_time=pcutils::hebing(otutab,metadata$Group)
#' time_by_cm(otu_time,n_cluster=4)->time_cm_res
#' plot.time_cm(time_cm_res)
#' }}
time_by_cm<-function(otu_time,n_cluster=6,min.std = 0){
  lib_ps("Mfuzz")

  otu_time=as.matrix(otu_time)
  mfuzz_class <- methods::new('ExpressionSet',exprs =otu_time)
  mfuzz_class <- Mfuzz::filter.NA(mfuzz_class, thres = 0.25)
  mfuzz_class <- Mfuzz::fill.NA(mfuzz_class, mode = 'mean')
  mfuzz_class <- Mfuzz::filter.std(mfuzz_class,min.std = min.std)
  #标准化数据
  mfuzz_class <- Mfuzz::standardise(mfuzz_class)
  set.seed(123)
  mfuzz_cluster <- Mfuzz::mfuzz(mfuzz_class, c = n_cluster, m = Mfuzz::mestimate(mfuzz_class))
  #mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = rownames(otu_time))

  cbind(mfuzz_class@assayData$exprs%>%as.data.frame(),
        id=mfuzz_class@assayData$exprs%>%as.data.frame()%>%rownames(),
        cluster=mfuzz_cluster$cluster,
        membership=mfuzz_cluster$membership%>%apply(.,1,max))->plotdat
  class(plotdat)=c("time_cm","data.frame")
  pcutils::del_ps("Mfuzz")
  return(plotdat)
}

#' Plot time_cm
#'
#' @param x time_cm
#' @param mem_thr membership threshold
#' @param ... add
#'
#' @return ggplot
#' @exportS3Method
#' @method plot time_cm
plot.time_cm=function(x,mem_thr=0.6,...){
  fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),
                  rep(0, length(c(255:150))))
  fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
  fancy.red <- c(c(0:255), rep(255, length(c(255:0))),
                 c(255:150))
  colo <- grDevices::rgb(b = fancy.blue/255, g = fancy.green/255,r = fancy.red/255)

  plotdat=x
  plotdat%>%filter(membership>=mem_thr)->plotdat1
  plotdat= reshape2::melt(plotdat1, id.vars = c("id","cluster","membership"),variable.name = "time")
  plotdat$cluster = factor(plotdat$cluster)

  ggplot(data=plotdat, aes(x=time, y=value, group=id,col=membership,alpha=membership)) +
    facet_wrap(cluster~.)+scale_color_gradientn(colours = colo)+
    geom_line(size=0.8)+ pctax_theme+
    scale_x_discrete(expand=c(0,0.2))+theme(plot.margin=unit(c(1,2,1,1),'lines'))
}


# mpse_da<-function(otutab,metadata_g,taxonomy,alpha=0.05){
#   lib_ps("MicrobiotaProcess")
#   data.frame(tax=taxonomy%>%apply(., 1, \(x)paste(unlist(x),collapse = "|")),otutab,check.names = F)->motu
#   write.table(motu,file = "./tmp",row.names = F,sep = "\t",quote = F)
#   if(ncol(metadata_g)!=2)stop("metadata_g need two columns, first is id, second is group")
#   colnames(metadata_g)[2]<-"Group"
#   MicrobiotaProcess::mp_import_metaphlan(profile="./tmp", mapfilename=metadata_g)->mpse
#   file.remove("./tmp")
#
#   mpse%>%
#     mp_diff_analysis(
#       .abundance = Abundance,
#       .group = Group,
#       first.test.alpha = alpha
#     )->mpse2
#
#   mpse2%<>%
#     mp_cal_abundance( # for each samples
#       .abundance = RareAbundance
#     ) %>%
#     mp_cal_abundance( # for each groups
#       .abundance=RareAbundance,
#       .group=Group
#     )
#
#   taxa.tree <- mpse2 %>%
#     mp_extract_tree(type="taxatree")
#   #taxa.tree %>% dplyr::select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Group, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))
#
#   p1<-mpse2%>%
#     mp_plot_diff_res(
#       group.abun = TRUE,
#       pwidth.abun=0.1,
#       tiplab.size=1
#     )
#
#   tibble::as_tibble(taxa.tree)%>%filter(!is.na(LDAmean))->saa1
#
#   tidyr::unnest(saa1,RareAbundanceBySample)->saa2
#
#   pp1<-ggplot(saa2, aes(x=label, y=RelRareAbundanceBySample,fill=Group)) +
#     geom_boxplot(aes(),width = 0.5) + ylab(label = "Abundance")+
#     #scale_fill_d3()+
#     coord_flip()+
#     ggfun::theme_stamp(
#       colour = c('grey90', 'white'),
#       axis = 'x',
#       axis.line.x = element_line(),
#       axis.title.y = element_blank(),legend.direction = "vertical")
#
#   pp2<-ggplot(saa1, aes(label, LDAmean)) +
#     ylab("LDA SCORE (log 10)")  +
#     geom_point(aes(col = Sign_Group))+
#     #scale_color_d3()+
#     #geom_bar(stat = "identity", aes(fill = Sign_Group),width = 0.5) + scale_fill_d3()+
#     coord_flip()+
#     ggfun::theme_stamp(
#       colour = c('grey90', 'white'),
#       axis = 'x',
#       axis.line.x = element_line(),
#       axis.title.y = element_blank(),
#       axis.line.y= element_blank(),
#       axis.ticks.y  = element_blank(),
#       axis.text.y=element_blank(),
#       legend.position = "none")
#   lib_ps("aplot",library = F)
#   p2<-pp1%>%insert_right(pp2,width = 0.6)
#
#   p3<-mpse2%>%
#     mp_plot_diff_cladogram(
#       label.size = 2.5,
#       hilight.alpha = .3,
#       bg.tree.size = .5,
#       bg.point.size = 2,
#       bg.point.stroke = .25
#     )
#
#   detach("package:MicrobiotaProcess")
#   return(list(tree=taxa.tree,p1=p1,p2=p2,p3=p3))
# }


#' Test the proper clusters k for c_means
#'
#' @param otu_group grouped otutab
#' @param filter_var filter the highest var
#' @param fast whether do the gap_stat?
#'
#' @return ggplot
#' @export
cm_test_k=function(otu_group,filter_var,fast=T){
  data_scaled=filter_top_var(otu_group,filter_var)

  #判断聚类个数
  #输入文件最好是按你想要的分组合并过的
  lib_ps("factoextra",library = F)
  #-------determining the number of clusters
  #1 Elbow method
  cp1<-factoextra::fviz_nbclust(data_scaled, kmeans, method = "wss",verbose = T) +
    labs(subtitle = "Elbow method")
  #2 Silhouette method
  cp2<-factoextra::fviz_nbclust(data_scaled, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")
  #3 Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  cp3=NULL
  if (!fast){
    cp3<-factoextra::fviz_nbclust(data_scaled, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
      labs(subtitle = "Gap statistic method")
  }
  return(list(cp1=cp1,cp2=cp2,cp3=cp3))
}

filter_top_var=function(otu_group,filter_var){
  #trans
  pcutils::dabiao("Filter top ",(1-filter_var)*100,  "% var and scale")
  group.var= apply(otu_group,1,var)
  otu_group.sel = otu_group[group.var >= quantile(group.var, filter_var),]#挑出变化较大的部分
  weight = c(apply(otu_group.sel, 1, var))
  data_scaled = pcutils::trans(otu_group.sel,method = 'standardize',MARGIN = 1)
  data_scaled
}

#' C-means cluster
#'
#' @param otu_group standardize data
#' @param k_num cluster number
#' @param filter_var filter the highest var
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab,package = "pcutils")
#' pcutils::hebing(otutab,metadata$Group)->otu_group
#' cm_test_k(otu_group,filter_var=0.7)
#' cm_res=c_means(otu_group,k=3,filter_var=0.7)
c_means<-function(otu_group,k_num,filter_var){
  lib_ps("e1071",library = F)

  data_scaled=filter_top_var(otu_group,filter_var)

  #-----Start clustering
  #set.seed(123)
  cm = e1071::cmeans(data_scaled, center=k_num, iter.max=500)

  cm_data = cbind.data.frame(Name = row.names(data_scaled),data_scaled,
                           Weight=apply(otu_group[rownames(data_scaled),], 1, var),
                           Cluster=cm$cluster,
                           Membership=apply(cm$membership, 1, max))
  res=list(data=otu_group,filter_var=filter_var,data_scaled=data_scaled,
           cm_data=cm_data,centers=cm$centers,membership=cm$membership)
  class(res)="cm_res"
  return(res)
}

#' Plot c_means result
#'
#' @param x a cm_res object
#' @param ... additional
#' @param filter_membership filter membership
#' @param mode 1~2
#' @param show.clust.cent show cluster center?
#' @param show_num show number of each cluster?
#'
#' @return ggplot
#' @exportS3Method
#' @method plot cm_res
plot.cm_res=function(x,filter_membership,mode=1,show.clust.cent=TRUE,show_num=TRUE,...){
  lib_ps("factoextra",library = F)

  cm_data=x$cm_data
  pcutils::dabiao('filter clusters, Membership >= ',filter_membership)
  cm_group = cm_data[cm_data$Membership>=filter_membership,]#筛选部分显著被聚类的项
  if(show_num){
    tmp=cm_group%>%count(Cluster)%>%mutate(new_cluster=paste0(Cluster,": ",n))
    cm_group$Cluster=setNames(tmp$new_cluster,tmp$Cluster)[cm_group$Cluster]
  }
  data_scaled=x$data_scaled

  if(mode==2){
    #show the cluster
    p<-factoextra::fviz_cluster(list(data = data_scaled[rownames(cm_group),], cluster=cm_group$Cluster),
                                   geom = c("point"),
                                   ellipse = TRUE,
                                   ellipse.alpha = 0.3, #used to be 0.6 if only points are plotted.
                                   ellipse.type = "norm",
                                   ellipse.level = 0.68,
                                   repel = TRUE,show.clust.cent = show.clust.cent
                                ) + pctax_theme
  }
  if(mode==1){

    cm_group.melt = reshape2::melt(cm_group, id.vars = c("Cluster","Membership", "Name","Weight"),variable.name = "Group")
    cm_group.melt$Cluster = factor(cm_group.melt$Cluster)

    p<-ggplot() +
      geom_line(data=cm_group.melt,
                aes(x=Group,y=value,group=Name,color=Cluster,alpha=Membership),
                linewidth=0.8)+
      pctax_theme +
      scale_x_discrete(expand=c(0,0))+theme(plot.margin=unit(c(1,2,1,1),'lines'))
    if(show.clust.cent){
      centers=x$centers%>%as.data.frame()

      if(show_num)centers$Cluster=tmp$new_cluster
      else centers$Cluster=rownames(centers)

      centers_dat=reshape2::melt(centers,id.vars = "Cluster",variable.name = "Group")

      p=p+
        scale_alpha_continuous(range = c(0.2,0.5))+
        geom_line(data = centers_dat,aes(x=Group,y=value,group=Cluster,color=Cluster),linewidth=3)
    }
  }
  cols1=pcutils::get_cols(length(unique(cm_group$Cluster)))
  return(p+scale_color_manual(values = cols1)+scale_fill_manual(values = cols1))
}
