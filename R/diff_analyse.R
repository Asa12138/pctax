#different analyse==============

bbtt<-function(pres_p,group){
  pres_p%>%arrange(p.format)%>%
    mutate(group=as.factor(Type),Tax=factor(Tax,levels = rev(Tax)))->pres1
  p1<-ggdotchart(pres1, x = "Tax", y = "p.format",color = "group",
                 dot.size = 7,  #palette = classcol, # 修改颜色
                 sorting = "ascending",
                 add = "segments",
                 add.params = list(color = "lightgray", size = 1.5),# 添加棒子
                 label = "Type",# 添加label
                 font.label = list(color = "white", size = 8,vjust = 0.5),
                 rotate = T,
                 ggtheme = theme_pubr()# 改变主题
  )+labs(x='')
  group=factor(group)
  pres1[,8:(7+nlevels(group))]%>%decostand(.,1,method = 'normalize')%>%
    mutate(Tax=rownames(.),Tax=factor(Tax,levels = rev(Tax)))->tmp
  p2<-  melt(tmp)%>%
    ggplot(.,aes(x=variable,y=Tax,fill=value))+geom_raster()+
    scale_fill_gradientn(colours = c("#1789E6", "#FFFFFF", "#B3192B"))+theme_pubr(legend = 'right')+
    theme_void()+
    theme(axis.text.x = element_text(),legend.position = 'top')
  p2+p1+plot_layout(widths = c(1.5,3))
}

ALDEX<-function(otutab,group){
  library(ALDEx2)
  if(nlevels(factor(group))==2){
    x.all <- aldex(otutab, group, mc.samples=16, test="t", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
  }
  if(nlevels(factor(group))>2){
    x.all <- aldex(otutab, group, mc.samples=16, test="kw", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
  }

  data.frame(Tax=rownames(x.all),x.all)->x.all
  x.all%>%mutate(p.signif=ifelse(glm.eBH>=0.05,'ns',ifelse(glm.eBH>=0.01,'*',ifelse(glm.eBH>0.001,'**','***'))))->x.all
  hebing(otutab,group)->tmp
  apply(tmp, 1, function(a)which(a==max(a))[[1]])->tmp$group
  apply(tmp, 1, function(a){return(colnames(tmp)[a[[4]]]) })->tmp$Type

  cbind(x.all,tmp[x.all$Tax,])->x.all
  cbind(x.all,otutab[x.all$Tax,])->x.all

  detach('package:ALDEx2')
  return(x.all)
}

ANCOM2<-function(otutab,metadata,group){
  source("D:/R/mycode/ANCOM-master/programs/ancom.R")

  res = ANCOM(otutab, metadata, struc_zero = NULL,main_var = group)

  res$out->out
  rownames(out)<-out$taxa_id
  hebing(otutab,metadata[,group])->tmp
  apply(tmp, 1, function(a)which(a==max(a))[[1]])->tmp$group
  apply(tmp, 1, function(a){return(colnames(tmp)[a[[4]]]) })->tmp$Type
  rename(out,Tax=taxa_id)->out
  cbind(out,tmp[out$Tax,])->out
  cbind(out,otutab[out$Tax,])->out

  # Step 3: Volcano Plot
  n_taxa = nrow(otutab)
  # Cutoff values for declaring differentially abundant taxa
  cut_off = c(0.9 * (n_taxa - 1), 0.8 * (n_taxa - 1), 0.7 * (n_taxa - 1), 0.6 * (n_taxa - 1))
  names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")
  # Annotation data
  dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
  res$fig = res$fig +
    geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
    geom_text(data = dat_ann, aes(x = x, y = y, label = label),
              size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
  return(res)
}

suijisenlin<-function(otutab,group,topN=10){
  suppressMessages(library(randomForest))
  t(otutab)%>%data.frame()%>%mutate(Group=group)%>%
    randomForest(Group ~.,data = .,importance=TRUE ,proximity=TRUE)->randomforest
  #plot(randomforest)
  #查看变量的重要性
  randomforest$importance%>%data.frame()%>%arrange(-MeanDecreaseAccuracy)->im
  im%>%head(topN)%>%mutate(Tax=rownames(.),Tax=factor(Tax,levels = rev(Tax)))->mim
  p1<-ggscatter(mim,x='MeanDecreaseAccuracy',y='Tax',color = '#7175E3',
                size="MeanDecreaseAccuracy",ggtheme = theme_pubr(legend = 'none'))+ylab('')

  hebing(otutab,group =group)->tmp
  tmp[as.character(mim$Tax),]%>%decostand(.,1,method = 'normalize')%>%
    mutate(Tax=rownames(.),Tax=factor(Tax,levels = rev(Tax)))->tmp
  p2<-melt(tmp)%>%
    ggplot(.,aes(x=variable,y=Tax,fill=value))+geom_raster()+
    scale_fill_gradientn(colours = c("#003366","white","#990033"))+theme_pubr(legend = 'right')+
    theme_void()+
    theme(axis.text.x = element_text(),legend.position = 'right')

  imp<-p1+p2+plot_layout(widths = c(3,1))

  return(list(im=im,imp=imp))
}

time_by_cm<-function(otu_time,n_cluster=6,mem_thr=0.6,min.std = 0){
  library(Mfuzz)
  fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),
                  rep(0, length(c(255:150))))
  fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
  fancy.red <- c(c(0:255), rep(255, length(c(255:0))),
                 c(255:150))
  colo <- rgb(b = fancy.blue/255, g = fancy.green/255,
              r = fancy.red/255)


  otu_time=as.matrix(otu_time)
  mfuzz_class <- new('ExpressionSet',exprs =otu_time)
  mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
  mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
  mfuzz_class <- filter.std(mfuzz_class,min.std = min.std)
  #标准化数据
  mfuzz_class <- standardise(mfuzz_class)
  set.seed(123)
  mfuzz_cluster <- mfuzz(mfuzz_class, c = n_cluster, m = mestimate(mfuzz_class))
  #mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = rownames(otu_time))

  cbind(mfuzz_class@assayData$exprs%>%as.data.frame(),
        id=mfuzz_class@assayData$exprs%>%as.data.frame()%>%rownames(),
        cluster=mfuzz_cluster$cluster,
        membership=mfuzz_cluster$membership%>%apply(.,1,max))->plotdat

  plotdat%>%filter(membership>=mem_thr)->plotdat1
  plotdat= melt(plotdat1, id.vars = c("id","cluster","membership"),
                variable.name = "time")
  plotdat$cluster = factor(plotdat$cluster)
  p=ggplot(data=plotdat, aes(x=time, y=value, group=id,col=membership,alpha=membership)) +
    facet_wrap(cluster~.)+scale_color_gradientn(colours = colo)+
    geom_line(size=0.8)+ theme_pubr(base_size = 12,legend = 'right')+
    scale_x_discrete(expand=c(0,0.2))+theme(plot.margin=unit(c(1,2,1,1),'lines'))
  return(list(time_res=plotdat1,p=p))
}

deseq_da<-function(otutab,metadata,group,ctrl=NULL){
  dplyr::select(metadata,!!group)%>%dplyr::rename(Group=1)->meta
  if(!is.null(ctrl))meta$Group<-relevel(meta$Group,ctrl)

  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = otutab, colData = meta, design = ~Group)#构建 DESeqDataSet 对象
  dds <- DESeq(dds) #差异分析

  #results(dds,name = resultsNames(dds)[2])%>%as.data.frame()
  res=data.frame()
  for (i in 2:length(resultsNames(dds))) {
    results(dds,name = resultsNames(dds)[i])%>%as.data.frame()->tmp
    tibble::rownames_to_column(tmp,"tax")->tmp
    res=rbind(res,data.frame(tmp,compare=resultsNames(dds)[i]))
  }
  res$compare=gsub("Group_","",res$compare)%>%gsub("_vs_","-",.)

  if(F){
    #limma
    #log后的数据哦
    list <- model.matrix(~factor(meta$Group)+0)  #把group设置成一个model matrix
    colnames(list) <- levels(factor(meta$Group))
    df.fit <- lmFit(otutab, list)  ## 数据与list进行匹配

    #df.matrix <- makeContrasts(KO - WT , levels = list)
    res=data.frame()
    for (i in 2:nlevels(meta$Group)){
      df.matrix <- makeContrasts(paste(levels(meta$Group)[c(i,1)],collapse = "-") , levels = list)
      fit <- contrasts.fit(df.fit, df.matrix)
      fit <- eBayes(fit)
      tempOutput <- topTable(fit,n = Inf)
      res=rbind(res,data.frame(tempOutput,compare=paste(levels(meta$Group)[c(i,1)],collapse = "-")))
    }
    colnames(res)=c("log2FoldChange","AveExpr","t","pvalue","padj","B","compare")
    rownames(res)->res$tax
  }

  logfc=0.5
  res[which(res$log2FoldChange >=logfc  & res$padj < 0.05),'sig'] <- 'up'
  res[which(res$log2FoldChange <= -logfc & res$padj < 0.05),'sig'] <- 'down'
  res[which(is.na(res$sig)),'sig'] <- 'none'
  res%>%mutate(tax1=ifelse(sig%in%c("up","down"),tax,""))->res
  res$label <- ifelse(res$sig%in%c("up","down"),"log2FoldChange>","log2FoldChange<")
  #unique(res$compare)
  #一对比较的火山图
  res%>%filter(!is.na(padj))->dat
  pp1<-ggplot(dat,aes(x=log2FoldChange,y=-log10(padj),color=sig))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=tax1),size=2)+
    scale_color_manual(values=c(up="#CC0000",none="#BBBBBB",down="#2f5688"),na.value ="#BBBBBB")+  #确定点的颜色
    facet_wrap(.~compare,scales = "free")+
    theme_bw()+  #修改图片背景
    theme(
      legend.title = element_blank()  #不显示图例标题
    )+
    ylab('-log10 (p-adj)')+  #修改y轴名称
    xlab('log2 (FoldChange)')+  #修改x轴名称
    geom_vline(xintercept=c(-logfc,logfc),lty=3,col="black",lwd=0.5) +  #添加垂直阈值|FoldChange|>2
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)  #添加水平阈值padj<0.05


  #多对比较的火山图
  res%>%filter(abs(log2FoldChange)>logfc)->dat
  res%>%group_by(compare)%>%summarise(max(log2FoldChange))%>%melt()->bardf
  res%>%group_by(compare)%>%summarise(min(log2FoldChange))%>%melt()->bardf1
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
                     width =0.4)+
    geom_text_repel(data=filter(dat,tax1!=""),aes(x = compare, y = log2FoldChange, label=tax1),
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
                       values = c("red","black"))

  pp2=p3+theme_minimal()+
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
  return(list(res=res,p1=pp1,p2=pp2))
}

mpse_da<-function(otutab,metadata_g,taxonomy,alpha=0.05){
  library(MicrobiotaProcess)
  data.frame(tax=taxonomy%>%apply(., 1, \(x)paste(unlist(x),collapse = "|")),otutab,check.names = F)->motu
  write.table(motu,file = "./tmp",row.names = F,sep = "\t",quote = F)
  if(ncol(metadata_g)!=2)stop("metadata_g need two columns, first is id, second is group")
  colnames(metadata_g)[2]<-"Group"
  mp_import_metaphlan(profile="./tmp", mapfilename=metadata_g)->mpse
  file.remove("./tmp")

  mpse%>%
    mp_diff_analysis(
      .abundance = Abundance,
      .group = Group,
      first.test.alpha = alpha
    )->mpse2

  mpse2%<>%
    mp_cal_abundance( # for each samples
      .abundance = RareAbundance
    ) %>%
    mp_cal_abundance( # for each groups
      .abundance=RareAbundance,
      .group=Group
    )

  taxa.tree <- mpse2 %>%
    mp_extract_tree(type="taxatree")
  #taxa.tree %>% dplyr::select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Group, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))

  p1<-mpse2%>%
    mp_plot_diff_res(
      group.abun = TRUE,
      pwidth.abun=0.1,
      tiplab.size=1
    )

  tibble::as_tibble(taxa.tree)%>%filter(!is.na(LDAmean))->saa1


  tidyr::unnest(saa1,RareAbundanceBySample)->saa2

  pp1<-ggplot(saa2, aes(x=label, y=RelRareAbundanceBySample,fill=Group)) +
    geom_boxplot(aes(),width = 0.5) + ylab(label = "Abundance")+
    #scale_fill_d3()+
    coord_flip()+
    ggfun::theme_stamp(
      colour = c('grey90', 'white'),
      axis = 'x',
      axis.line.x = element_line(),
      axis.title.y = element_blank(),legend.direction = "vertical")

  pp2<-ggplot(saa1, aes(label, LDAmean)) +
    ylab("LDA SCORE (log 10)")  +
    geom_point(aes(col = Sign_Group))+
    #scale_color_d3()+
    #geom_bar(stat = "identity", aes(fill = Sign_Group),width = 0.5) + scale_fill_d3()+
    coord_flip()+
    ggfun::theme_stamp(
      colour = c('grey90', 'white'),
      axis = 'x',
      axis.line.x = element_line(),
      axis.title.y = element_blank(),
      axis.line.y= element_blank(),
      axis.ticks.y  = element_blank(),
      axis.text.y=element_blank(),
      legend.position = "none")

  p2<-pp1%>%insert_right(pp2,width = 0.6)

  # p2<-mpse2%>%
  #   mp_plot_diff_boxplot(
  #     .group = Group
  #   )

  p3<-mpse2%>%
    mp_plot_diff_cladogram(
      label.size = 2.5,
      hilight.alpha = .3,
      bg.tree.size = .5,
      bg.point.size = 2,
      bg.point.stroke = .25
    )

  detach("package:MicrobiotaProcess")
  return(list(tree=taxa.tree,p1=p1,p2=p2,p3=p3))
}

cm_test_k<-function(otutab,group,fast=T){
  #1数据处理
  hebing(otutab,group)->date_gen_cpm
  season.var = apply(date_gen_cpm,1,var)
  var.thre = quantile(season.var, 0.75)#挑出变化较大的部分，阈值
  df.season.sel = date_gen_cpm[season.var >= var.thre,]
  weight = c(apply(df.season.sel, 1, var))
  wtf = decostand(df.season.sel,method = 'standardize',MARGIN = 1)

  #2判断聚类个数
  #输入文件最好是按你想要的分组合并过的
  library(cluster)
  library(factoextra)
  #-------determining the number of clusters
  #1 Elbow method
  cp1<-fviz_nbclust(wtf, kmeans, method = "wss") +
    labs(subtitle = "Elbow method")
  #2 Silhouette method
  cp2<-fviz_nbclust(wtf, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")
  # Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  #3 Use verbose = FALSE to hide computing progression.
  #set.seed(123)
  cp3=NULL
  if (!fast){
    cp3<-fviz_nbclust(wtf, kmeans, nstart = 25,  method = "gap_stat", nboot = 100)+
      labs(subtitle = "Gap statistic method")
  }
  return(list(cp1=cp1,cp2=cp2,cp3=cp3,wtf=wtf,weight=weight))
}

c_means<-function(otutab,wtf,k,weight){
  library(e1071)
  library(NbClust)
  #-----Start clustering
  #set.seed(123)
  cm = cmeans(wtf, center=k, iter.max=500)
  #cm$cluster = factor(cm$cluster, levels=c(1,2,3,4))
  cat('被聚类的项')
  table(cm$cluster)%>%print()
  #show the cluster
  cmp1<-fviz_cluster(list(data = otutab[rownames(wtf),], cluster=cm$cluster),
                     geom = c("point"),
                     ellipse = TRUE,
                     ellipse.alpha = 0.3, #used to be 0.6 if only points are plotted.
                     ellipse.type = "norm",
                     ellipse.level = 0.68,
                     repel = TRUE) + theme_pubr(legend = 'right')

  tempp = cbind.data.frame(wtf, Weight=weight, Cluster=cm$cluster, Membership=apply(cm$membership, 1, max), Taxon = row.names(wtf))
  cm_group = tempp[tempp$Membership>=0.65,]#筛选部分显著被聚类的项
  cat('筛选部分显著被聚类的项,Membership>=0.65')
  table(cm_group$Cluster)%>%print()
  cm_group.melt = melt(cm_group, id.vars = c("Cluster","Membership", "Taxon","Weight"),
                       variable.name = "Date")
  cm_group.melt$Cluster = factor(cm_group.melt$Cluster)

  cmp2<-ggplot(data=cm_group.melt, aes(x=Date, y=value, group=Taxon, color=Cluster, size=Weight,alpha=Membership)) +
    geom_line(size=0.8)+ theme_pubr(base_size = 12,legend = 'right') +
    scale_x_discrete(expand=c(0,0))+theme(plot.margin=unit(c(1,2,1,1),'lines'))+
    scale_color_npg()
  #scale_y_continuous(limits = c(-1.5,2))+
  #theme(legend.position = "none") #expand c(0,0) eliminates the sapcing4

  return(list(cmp1=cmp1,cmp2=cmp2,cm_group=cm_group))
}

multibox<-function(clusterpvalue,tax,group){
  pols=list()
  for (i in tax){
    Tax=i
    #使用log做scale，可以减少离群值的影响
    plotdat<-data.frame(value=as.vector(decostand(t(clusterpvalue[Tax,(10+nlevels(group)):(ncol(clusterpvalue)-2)]),method = 'log',MARGIN = 2)),
                        Group=factor(group))

    if(F){
      plotdat$Group.1<-as.numeric(plotdat$Group)
      #错在 group 的数据类型，他是因子型，但如果想将坐标点连接起来得为数值
      aggregate(plotdat$value,by=list(plotdat$Group1),median)->a
      #加line的箱型图
      ggplot(data = plotdat,aes(x=Group.1,y=value))+
        geom_boxplot(aes(group=Group,color=Group),outlier.shape = NA)+
        geom_line(data = a,aes(y=x),color='skyblue',size=1.2)+
        #geom_smooth(size=1)+
        #stat_compare_means(show.legend = FALSE)+
        scale_x_continuous(breaks=unique(plotdat$Group.1),labels = levels(plotdat$Group))+
        labs(y=Tax,x=NULL)+
        geom_jitter(aes(group=Group,color=Group),width = 0.15,alpha=0.8,size=1)+
        scale_color_npg(name='Date')+
        #scale_color_manual(values =brewer.pal(12,'Paired'))+
        theme(plot.margin=unit(rep(0.5,4),'lines'))+
        theme_pubr(legend = 'right',x.text.angle = 45)
    }
    pols[[i]]=ggplot(data = plotdat,aes(x=Group,y=value))+
      geom_boxplot(aes(group=Group,color=Group),outlier.shape = NA)+
      #geom_line(data = a,aes(y=x),color='skyblue')+
      #geom_smooth(size=1)+
      stat_compare_means(show.legend = FALSE,label.x = 1)+
      labs(y=Tax,x=NULL)+
      geom_jitter(aes(group=Group,color=Group),width = 0.15,alpha=0.8,size=1)+
      #scale_color_manual(name='Season',values = seacol)+
      scale_color_manual(values =brewer.pal(12,'Paired'))+
      theme(plot.margin=unit(rep(0.5,4),'lines'))+
      theme_pubr(legend = 'none')

  }
  return(pols)

}
