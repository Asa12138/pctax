#b_diversity============

#' Lm for sample similarity and geographical distance
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param geo a two-columns dataframe, first is latitude, second is longitude
#' @param ... additional
#' @import NST
#' @return a ggplot
#' @export
#' @references 1. Graco-Roza, C. et al. Distance decay 2.0 - A global synthesis of taxonomic and functional turnover in ecological communities. Glob Ecol Biogeogr 31, 1399–1421 (2022).
#' @examples
#'data(otutab)
#'data.frame(row.names = rownames(metadata),lat=runif(18,30,35),long=runif(18,40,45))->geo
#'geo_sim(otutab,geo)
geo_sim<-function(otutab,geo,method="bray",spe_nwk=NULL,...){
  lib_ps("NST","geosphere")
  #经纬度数据转换
  #直接欧式距离算不太准
  #pcutils::toXY(geo)%>%ggscatter(.,"X","Y",label = rownames(geo))
  geosphere::distm(geo[,2:1])->geo_dist
  rownames(geo_dist)<-colnames(geo_dist)<-rownames(geo)

  geo_dist%>%as.dist()%>%NST::dist.3col()%>%mutate(dis=dis/1000)->geo_dist
  #这一步可以选择各种b_diversity指数，甚至是谱系与功能多样性指数
  b_dist(otutab,method,spe_nwk)%>%NST::dist.3col()%>%mutate(dis=1-dis)->similarity

  merge(geo_dist,similarity,by = c("name1","name2"),suffixes = c(".geo",".b"))->a
  pcutils::my_lm(a[,4:3],...)+xlab("Distance(km)")+ylab(paste0("1- ",method))
}

#' Beta_diversity Ordination: dimensionality reduction
#'
#' @param object object
#' @param ... additional
#'
#' @return b_res object
#' @export
#' @references \link{https://www.jianshu.com/p/9694c0b6302d}
#'\link{https://mixomicsteam.github.io/Bookdown/plsda.html}
#'\link{https://zhuanlan.zhihu.com/p/25501130}
#' @examples
#' data(otutab)
#' b_analyse(otutab,method = "pca")->b_res
#' plot(b_res,metadata$Group)
#' #b_analyse(otutab,method = "all",group=metadata$Group)->b_res
#' #plot(b_res,metadata$Group,mode=3)
b_analyse <- function(object,...){
  UseMethod("b_analyse", object)
}

#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param norm should normalized or not? (hellinger)
#' @param method one of "pca","pcoa","ca","dca","nmds","plsda","tsne","umap","lda","all"
#' @param group if needed, give a group vector
#' @param dist if use pcoa, your can choose a dist method(default: euclidean)
#' @param ndim how many dimension be kept?(default:2) 3 for b_res_3d()
#'
#' @exportS3Method
#' @rdname b_analyse
b_analyse.data.frame<-function(otutab,norm=T,method=c("pca","nmds"),group=NULL,dist = 'euclidean',ndim=2){
  lib_ps("ade4","vegan","dplyr")
  all=c("pca","pcoa","ca","dca","nmds","plsda","tsne","umap","lda","all")
  if(!all(method%in%all))stop("method is one of ",paste(all,collapse = ","))
  if("all"%in%method)method=all[-length(all)]
  data.frame(t(otutab))->dat
  if(norm)dat.h <- decostand(dat, "hellinger") else dat.h<-dat
  #storage dataframes
  data.frame(name=colnames(otutab))->sample_site
  data.frame(eigs=paste('eig',1:ndim))->sample_eig
  data.frame(var=rownames(otutab))->var_site
  data.frame(var=rownames(otutab))->var_contri
  #PCA:principal components analysis 主成分分析
  if("pca"%in%method){
    # dat.pca <- rda(dat.h, scale = F)#scale=T,对列变量进行标准化
    # summary(dat.pca, scaling = 1)->a#scaling=1展示样本距离，=2展示变量夹角
    # cbind(sample_eig,a[["cont"]][["importance"]][2,1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'PCA'
    # cbind(sample_site,a$sites[,1:2])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('PCA1', 'PCA2')

    dat.pca = ade4::dudi.pca(dat.h, scale=F, scannf=FALSE, nf=3)
    cbind(sample_eig,(dat.pca$eig/sum(dat.pca$eig))[1:ndim])->sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- 'PCA'
    cbind(sample_site,dat.pca$li[,1:ndim])->sample_site#绘图坐标
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("PCA",1:ndim)
    cbind(var_site,dat.pca$c1[,1:2])->var_site#变量坐标
    colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- paste0("PCA",1:2)
    dat.pca$co[,1:2]%>%
      transmute(contri=((Comp1^2/(dat.pca$eig)[1])+(Comp2^2/(dat.pca$eig)[2]))*100/2)%>%
      cbind(var_contri,.)->var_contri
    colnames(var_contri)[ncol(var_contri)] <- 'PCA'

    #这个结果可以用fviz_pca直接可视化
    if(F){
      fviz_screeplot(dat.pca)
      fviz_pca_ind(dat.pca, col.ind = "cos2")
      fviz_pca_var(dat.pca)
    }
  }
  #PcoA:principal coordinate analysis 主坐标分析
  if("pcoa"%in%method){
    # dat.pcoa <- ape::pcoa(vegdist(dat.h,method = dist))
    # cbind(sample_eig,dat.pcoa$values$Relative_eig[1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'PCoA'
    # cbind(sample_site,dat.pcoa$vectors[,c(1,2)])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('PCoA1', 'PCoA2')

    dat.pco= ade4::dudi.pco(vegan::vegdist(dat.h,method = dist), scannf=FALSE, nf=3)
    cbind(sample_eig,(dat.pco$eig/sum(dat.pco$eig))[1:ndim])->sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- 'PCoA'
    cbind(sample_site,dat.pco$li[,1:ndim])->sample_site#绘图坐标
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("PCoA",1:ndim)
    # cbind(var_site,dat.pco$c1[,1:2])->var_site#变量坐标
    # colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('PCoA1', 'PCoA2')
    # dat.pco$co[,1:2]%>%transmute(contri=((Comp1^2/(dat.pco$eig)[1])+(Comp2^2/(dat.pco$eig)[2]))*100/2)%>%
    #   cbind(var_contri,.)->var_contri
    # colnames(var_contri)[ncol(var_contri)] <- 'PCoA'
  }
  #CA:correspondence analysis 对应分析
  if("ca"%in%method){
    # dat.ca <- cca(dat.h,scale = F)
    # summary(dat.ca, scaling = 1)->b
    # cbind(sample_eig,(dat.ca$CA$eig/sum(dat.ca$CA$eig))[1:2])->sample_eig#两轴的解释度
    # colnames(sample_eig)[ncol(sample_eig)] <- 'CA'
    # cbind(sample_site,b$sites[,1:2])->sample_site#绘图坐标
    # colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('CA1', 'CA2')

    dat.coa= ade4::dudi.coa(dat.h, scannf=FALSE, nf=3)
    cbind(sample_eig,(dat.coa$eig/sum(dat.coa$eig))[1:ndim])->sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- 'CA'
    cbind(sample_site,dat.coa$li[,1:ndim])->sample_site#绘图坐标
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("CA",1:ndim)
    cbind(var_site,dat.coa$co[,1:2])->var_site#变量坐标
    colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('CA1', 'CA2')
    dat.coa$c1[,1:2]%>%transmute(contri=((CS1^2/(dat.coa$eig)[1])+(CS2^2/(dat.coa$eig)[2]))*100/2)%>%
      cbind(var_contri,.)->var_contri
    colnames(var_contri)[ncol(var_contri)] <- 'CA'

  }
  #DCA:detrended correspondence analysis 去趋势分析
  if("dca"%in%method){
    dat.dca= vegan::decorana(dat.h)
    #all.dca=summary(dat.dca)
    cbind(sample_eig,(dat.dca$evals/sum(dat.dca$evals))[1:ndim])->sample_eig
    colnames(sample_eig)[ncol(sample_eig)] <- 'DCA'
    cbind(sample_site,dat.dca$rproj[,1:ndim])->sample_site#绘图坐标
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("DCA",1:ndim)
    # cbind(var_site,dat.pco$c1[,1:2])->var_site#变量坐标
    # colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('DCA1', 'DCA2')
    # dat.pco$co[,1:2]%>%transmute(contri=((Comp1^2/(dat.pco$eig)[1])+(Comp2^2/(dat.pco$eig)[2]))*100/2)%>%
    #   cbind(var_contri,.)->var_contri
    # colnames(var_contri)[ncol(var_contri)] <- 'DCA'
  }
  #NMDS:non-metric multi-dimensinal scaling 非度量多维尺度分析
  if("nmds"%in%method){
    suppressMessages(dat.nmds <- vegan::metaMDS(dat.h,k = ndim, distance = dist))
    cbind(sample_eig,data.frame(a=rep("",ndim)))->sample_eig#两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- 'NMDS'
    #cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site,dat.nmds$points)->sample_site
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("NMDS",1:ndim)
    cbind(var_site,dat.nmds$species[,1:2])->var_site#变量坐标
    colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('NMDS1', 'NMDS2')
    cbind(var_contri,var_site%>%transmute(((NMDS1^2+NMDS2^2)/2)))->var_contri
    colnames(var_contri)[ncol(var_contri)] <- 'NMDS'

  }
  #PLS-DA:partial least squares discriminant analysis 偏最小二乘法判别分析
  if("plsda"%in%method){
    if(is.null(group))stop("plsda need group!")
    lib_ps("mixOmics")
    mixOmics::plsda(dat.h,group,ncomp = ndim)->plsdat
    cbind(sample_eig,data.frame(plsdat=plsdat$prop_expl_var$X))->sample_eig#两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- 'PLS_DA'
    #cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site,plsdat$variates$X)->sample_site
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("PLS_DA",1:ndim)
    cbind(var_site,plsdat$mat.c[,1:2]*100)->var_site#变量坐标
    colnames(var_site)[(ncol(var_site)-1):ncol(var_site)] <- c('PLS_DA1', 'PLS_DA2')
    var_site[,(ncol(var_site)-1):ncol(var_site)]%>%transmute(contri=((PLS_DA1^2/(plsdat$prop_expl_var$X)[1])+(PLS_DA2^2/(plsdat$prop_expl_var$X)[2]))*100/2)%>%
      cbind(var_contri,.)->var_contri
    colnames(var_contri)[ncol(var_contri)] <- 'PLS_DA'
  }
  #LDA:linear discriminant analysis 线性判别分析
  if("lda"%in%method){
    if(is.null(group))stop("lda need group!")
    lda.sol = MASS::lda(dat.h, group)
    P = lda.sol$scaling
    # 将均值向量降维
    means = lda.sol$means %*% P
    # 加权平均的出总的降维均值向量，权重就是lda.sol$prior
    total_means = as.vector(lda.sol$prior %*% means)
    # 把样本降维并平移
    x<-as.matrix(dat.h) %*% P - (rep(1, nrow(dat.h)) %o% total_means)
    #ldat=(lda.sol$svd)**2/sum((lda.sol$svd)**2)#两轴的解释度
    cbind(sample_eig,a="")->sample_eig#两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- 'LDA'
    #cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site,x)->sample_site
    colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('LDA1', 'LDA2')
  }
  #t-SNE:t-distributed stochastic neighbor embedding t-分布式随机邻域嵌入
  if("tsne"%in%method){
    lib_ps("Rtsne")
    Rtsne::Rtsne(dat.h,dims = ndim,pca = T,
                 max_iter = 1000,
                 theta = 0.4,
                 perplexity = 5,#默认20，会报错
                 verbose = F)->tsne
    cbind(sample_eig,a="")->sample_eig#两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- 't_SNE'
    #cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site,tsne$Y)->sample_site
    colnames(sample_site)[(ncol(sample_site)-ndim+1):ncol(sample_site)] <- paste0("t_SNE",1:ndim)
  }
  #UMAP:uniform manifold approximation and projection 均一流形近似投影
  if("umap"%in%method){
    lib_ps("umap")
    umap::umap(dat.h)->umapr
    cbind(sample_eig,a="")->sample_eig#两轴的解释度
    colnames(sample_eig)[ncol(sample_eig)] <- 'UMAP'
    #cbind(sample_site,scores(dat.nmds))->sample_site
    cbind(sample_site,umapr$layout)->sample_site
    colnames(sample_site)[(ncol(sample_site)-1):ncol(sample_site)] <- c('UMAP1', 'UMAP2')
  }
  #fso: fuzzy set ordination 模糊集排序
  if(F){
    lib_ps("fso")
    env.fso<-fso(metadata[,12:13],vegdist(dat.h))
    data.frame(x=env.fso$data,y=env.fso$mu)%>%ggplot(.,aes(x=y.1,y=y.2))+geom_point(aes(col=metadata$Group))
  }
  #sofm:self-organizing feature map 自组织特质
  if(F){
    lib_ps("kohonen")
    kohonen::xyf()
  }

  print('four dataframes in a list, 1 is eig, 2 is sample_site, 3 is var, 4is var contribution')
  b_res=list(sample_eig=sample_eig,sample_site=sample_site,var_site=var_site,var_contri=var_contri)
  class(b_res)="b_res"
  return(b_res)
}


#' @exportS3Method
#' @rdname b_analyse
b_analyse.pc_otu<-function(pc,tbl="otutab",norm=T,method=c("pca","ca"),group=NULL,dist = 'euclidean'){
  pc_valid(pc)
  otutab<-pc$tbls[[tbl]]
  pc$b_res<-b_analyse.data.frame(otutab,norm=norm,method=method,group=group,dist = dist)
  return(pc)
}

#' base plot for pca/rda
#'
#'
plot_b_like<-function(plotdat,mode=1,pal=NULL,sample_label=T){
  if (mode==1){
    plist <- {ggplot(plotdat, aes(x=x1, y=x2))+
        geom_point(aes(bg=level),pch = 21, colour = "black",size = 2)+ #可在这里修改点的透明度、大小
        geom_vline(xintercept = 0, color = 'gray', linewidth = 0.4) +
        geom_hline(yintercept = 0, color = 'gray', linewidth = 0.4)}
    if(!is.numeric(plotdat$level)){plist=plist+
      stat_ellipse(aes(fill=level),type="norm",geom="polygon",alpha=0.1,color=NA,level = 0.68)}
  }
  else if (mode==2){
    plist <- {ggplot(plotdat, aes(x=x1, y=x2))+
        geom_point(aes(color=level),size = 2)+ #可在这里修改点的透明度、大小
        geom_vline(xintercept = 0, color = 'gray', linewidth = 0.4) +
        geom_hline(yintercept = 0, color = 'gray', linewidth = 0.4)}
    if(!is.numeric(plotdat$level))plist=plist+stat_ellipse(aes(color=level),level = 0.68)
  }
  else if (mode==3){
    if(is.numeric(plotdat$level))stop("Group is continous!")
    centroid <- aggregate(cbind(x1,x2) ~ level,
                          data = plotdat,FUN = mean)
    plotdat1<-dplyr::left_join(plotdat, centroid, by = "level",suffix = c("",".cen"))
    plotdat<-data.frame(plotdat1,row.names = rownames(plotdat))
    plist <- {ggplot(plotdat, aes(x=x1, y=x2))+
        geom_point(aes(color=level),size = 2)+ #可在这里修改点的透明度、大小
        geom_vline(xintercept = 0, color = 'gray', linewidth = 0.4) +
        geom_hline(yintercept = 0, color = 'gray', linewidth = 0.4) +
        geom_segment(aes(xend=x1.cen,yend=x2.cen,color=level),show.legend = F)+
        geom_label(data = centroid, aes(label = level, fill = level), size = 5,
                   show.legend = FALSE,color="white")
    }
  }
  #sample_label
  if (sample_label){
    lib_ps("ggrepel")
    plist=plist+geom_text_repel(aes(x=x1, y=x2,label=rownames(plotdat)),
                                          col='black',size=2.5)
  }

  if(!is.numeric(plotdat$level)){plist=plist+
    scale_color_manual(values = pal) + #可在这里修改点的颜色
    scale_fill_manual(values = pal)}
  else {plist=plist+
    scale_color_gradientn(colours = brewer.pal(8,"Reds"))
  }
  plist=plist+theme_classic()
  return(plist)
}

#' Plot a b_res
#'
#' @param b_res a b_res object
#' @param Group group vector for color
#' @param mode plot mode:1~3
#' @param bi plot variables segments?
#' @param Topn how many variables to show?
#' @param rate segments length rate
#' @param margin plot the margin boxplot?
#' @param box margin plot box or density?
#' @param pal colors for group
#' @param sample_label plot the labels of samples?
#'
#' @return a ggplot
#' @exportS3Method
#'
#' @seealso \code{\link{b_analyse}}
#'
plot.b_res<-function(b_res,Group,mode=1,bi=F,Topn=10,rate=1,margin=F,box=T,pal=NULL,sample_label=T){
  lib_ps("dplyr","ggplot2","ggnewscale")
  if(is.null(pal)&!is.numeric(Group))pal=pcutils::get_cols(n = length(unique(Group)),pal = RColorBrewer::brewer.pal(5,"Set2"))
  #mode 代表用哪种风格画图，常用的1-3已经准备好了，临时改的话添加4就行。
  plist=list()
  for (i in colnames(b_res$sample_eig)[-1]){
    b_res$sample_site%>%dplyr::select(starts_with(i))->tmp
    plotdat<-data.frame(x1=tmp[,1],x2=tmp[,2],level=Group,row.names = b_res$sample_site$name)
    b_res$sample_eig%>%dplyr::select(starts_with(i))%>%unlist->eig

    plist[[i]]=plot_b_like(plotdat,mode=mode,pal=pal,sample_label=sample_label)
    #labs on axis
    if(i%in%c("PCA","PCoA","CA","PLS_DA")){
      plist[[i]] <- plist[[i]]+
        labs(x = paste(paste(i,'1: ',sep = ''), round(100 * eig[1], 2), '%'),
             y = paste(paste(i,'2: ',sep = ''), round(100 * eig[2], 2), '%'))
    }
    else {
      plist[[i]] <- plist[[i]]+labs(x = paste0(i,"1"),y = paste0(i,"2"))
    }

    if(bi==T){
      b_res$var_site%>%dplyr::select(starts_with(i))->tmp
      tmp*rate->tmp
      if((b_res$var_site%>%dplyr::select(starts_with(i))%>%ncol())>0){
        cbind(var=b_res$var_site[,1],tmp,contri=b_res$var_contri[,i])->bi_tmp
        colnames(bi_tmp)<-c("var","x1","x2","contri")
        bi_tmp%>%top_n(Topn,contri)->bi_tmp1
        plist[[i]] <- plist[[i]]+
          ggnewscale::new_scale_color()+
          geom_segment(data = bi_tmp1, aes(x = 0, y = 0, xend = x1,yend = x2,color=contri),
                       size = 0.3,arrow = arrow(length=unit(0.1, "inches"))) +
          geom_text(data = bi_tmp1, aes(x =  x1*1.05, y = x2*1.05,
                                        color=contri, label = var), size = 2)+
          scale_color_gradientn(name="contribution",colours = c("#00AFBB", "#E7B800", "#FC4E07"))
      }
    }

    if(margin){
      lib_ps("aplot")
      if(box){
      p1<-ggplot(plotdat,aes(x=level,y=x2,fill=level))+
        geom_boxplot(outlier.shape = NA)+
        ylab(label = NULL)+xlab(label = NULL)+
        scale_fill_manual(values = pal)+
        theme_transparent(base_size = 11)+
        theme(legend.position = "none",
              axis.line.x = element_line(),
              axis.ticks.x = element_line(),
              axis.text.x = element_text(angle = 90,vjust = 0.5,size = rel(0.8)))

      p2<-ggplot(plotdat,aes(x=level,y=x1,fill=level))+
        geom_boxplot(outlier.shape = NA)+
        ylab(label = NULL)+xlab(label = NULL)+
        scale_fill_manual(values = pal)+
        theme_transparent(base_size = 11)+
        theme(legend.position = "none",
              axis.line.y = element_line(),
              axis.ticks.y = element_line(),
              axis.text.y = element_text(size = rel(0.8)))+coord_flip()
      }
      else{
      p1<-ggplot(plotdat,aes(y=x2,fill=level,col=level))+
        geom_density(alpha=0.5)+
        ylab(label = NULL)+xlab(label = NULL)+
        scale_fill_manual(values = pal)+
        scale_color_manual(values = pal)+
        theme_transparent()+theme(legend.position = "none")
      p2<-ggplot(plotdat,aes(x=x1,fill=level,col=level))+
        geom_density(alpha=0.5)+
        ylab(label = NULL)+xlab(label = NULL)+
        scale_fill_manual(values = pal)+
        scale_color_manual(values = pal)+
        theme_transparent()+theme(legend.position = "none")
      }
      plist[[i]]%>%insert_right(p1,width = 0.2)%>%insert_top(p2,height = 0.2)->plist[[i]]
    }

  }
  return(plist)
}


#' 3D plot for b_res
#'
#' @param b_res a b_res object
#' @param group group vector for color
#'
#' @return plotly list
#' @export
#'
#' @examples
#' data(otutab)
#' b_analyse(otutab,method = "pca",ndim=3)->b_res
#' b_res_3d(b_res,metadata$Group)
b_res_3d<-function(b_res,group){
  lib_ps("plotly")
  plist=list()
  for (i in colnames(b_res$sample_eig)[-1]){
    b_res$sample_site%>%dplyr::select(starts_with(i))->tmp
    if(ncol(tmp)!=3)next
    plotdat<-data.frame(dim1=tmp[,1],dim2=tmp[,2],dim3=tmp[,3],level=factor(group))
    plot_ly(plotdat,x=~dim1,y=~dim2,z=~dim3,color = ~level,type ="scatter3d", mode ="markers")%>%
      plotly::layout(title=i)->plist[[i]]
  }
  return(plist)
}


#' Permanova between a otutab and a variable
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param envs factors need to test
#' @param norm should normalize?(default:T)
#' @param method adonis/mrpp/anosim/mantel
#' @param two two by two adonis test
#' @param each test factor one by one, rather than whole
#'
#' @return a g_test object with these columns
#' \item{group}{the test group or factor}
#' \item{r}{relationship}
#' \item{r2}{model R-square}
#' \item{p_value}{model test p_value}
#' \item{sig}{whether significant}
#' @export
#' @references \link{https://blog.csdn.net/qq_42458954/article/details/110390488}
#' @examples
#' data(otutab)
#' permanova(otutab,metadata[,c(2:4,12:14),drop=F])->adonis_res
#' sanxian(adonis_res)
#' plot(adonis_res)
permanova<-function(otutab,envs,norm=T,each=T,method="adonis",two=F){
  lib_ps("vegan")
  all=c("adonis","anosim","mrpp","mantel")
  if(!method%in%all)stop(paste0("method should be one of ",all))
  set.seed(123)
  env=envs
  data.frame(t(otutab))->dat
  if(norm)otu.t <- decostand(dat, "hellinger") else otu.t<-dat
  if(each){
    #adnois不只是检验分组变量，连续的数值变量也可以检验。
    soil <- NULL
    #print(paste0('==========',colnames(env)[i],' adnois result=========='))
    if(method=="adonis"){
      for (i in 1:ncol(env)){
        dat.div <- adonis2(otu.t ~(env[,i]),permutations = 999, method="bray")
        soil <- rbind(soil,c(colnames(env)[i],dat.div$R2[1], dat.div$`Pr(>F)`[1]))
        if(two){if((is.factor(env[,i])|class(env[,i])=="Date"|is.character(env[,i]))){
          env[,i]%>%as.factor()->group
          lib_ps(pairwiseAdonis)
          dat.pairwise.adonis <- pairwise.adonis(x=otu.t, factors=group, sim.function = "vegdist",
                                                 sim.method = "bray",p.adjust.m = "BH",
                                                 reduce = NULL,perm = 999)
          sanxian(dat.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")],rows=NULL)%>%print()
        }}
      }
    }
    if(method=="mantel"){
      #only numeric variables can do with mantel
      env%>%select_if(\(x)is.numeric(x)&!is.factor(x))->env
      species.distance<-vegdist(otu.t,method = 'bray')
      for (i in 1:ncol(env)){
        dd <- mantel(species.distance, vegdist(env[,i], method = "bray"),
                     method = "pearson", permutations = 999, na.rm = TRUE)
        soil <- rbind(soil,c(colnames(env)[i],dd$statistic, dd$signif))
      }
    }
    #only group variables can do with mrpp/anosim
    if(method=="anosim"){
      env%>%select_if(\(x)class(x)=="Date"|is.factor(x)|is.character(x))->env
      env%>%mutate_all(\(x)as.factor(x))->env
      for (i in 1:ncol(env)){
        env[,i]->group
        anosim_res=anosim(otu.t,group,permutations = 999)
        soil <- rbind(soil,c(colnames(env)[i],anosim_res$statistic, anosim_res$signif))
      }
    }
    if(method=="mrpp"){
      env%>%select_if(\(x)class(x)=="Date"|is.factor(x)|is.character(x))->env
      env%>%mutate_all(\(x)as.factor(x))->env
      for (i in 1:ncol(env)){
        env[,i]->group
        mrpp_res=mrpp(otu.t,group,permutations = 999)
        soil <- rbind(soil,c(colnames(env)[i],mrpp_res$A, mrpp_res$Pvalue))
      }
    }
    soil <- data.frame(group=soil[,1],r2=soil[,2],p_value=soil[,3])
  }
  else{
    dat.div <- adonis2(otu.t ~.,data = env,permutations = 999, method="bray")
    dat.div<-dat.div[seq_len(nrow(dat.div)-2),,drop=F]%>%as.data.frame()
    soil<-data.frame(group=rownames(dat.div),r2=dat.div$R2,p_value=dat.div$`Pr(>F)`)
  }
  soil$r2=round(as.numeric(soil$r2),4)
  soil$p_value=round(as.numeric(soil$p_value),4)
  soil$sig=soil$p_value <0.05
  if(method%in%c("anosim","mrpp","mantel"))colnames(soil)[2]="r"
  class(soil)<-c("g_test","data.frame")
  return(soil)
}

#' Plot g_test
#'
#' @param aa a g_test object
#'
#' @return ggplot
#' @exportS3Method
#'
#' @seealso  \code{\link{permanova}}
plot.g_test<-function(aa){
  if("r2"%in%colnames(aa)){
    aa$group = factor(aa$group,levels = aa$group[order(aa$r2)])#按R2值排序
    aa$sig_l=cut(aa$p_value, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("**", "*", ""))
    p<-ggplot(aa, aes(x =group, y = r2),size=2) +
      geom_bar(stat = 'identity', width = 0.8,color="black",fill="#5AAFD1") +
      scale_fill_manual(guide = 'none')+labs(title = "Environmental factor")+
      geom_text(aes(y = r2*0.9, label = sig_l),size = 5) +#可调星号位置
      xlab(NULL)+ylab(expression(r^"2"))+
      scale_y_continuous(expand = c(0,0))+
      theme_classic()
  }
  if("r"%in%colnames(aa)){
    aa$sig_l=ifelse(aa$sig,"sig","nosig")
    p<-ggplot(data=aa,aes(x=group,y=r),size=2)+
      geom_point(aes(size=r,color=sig_l),alpha=0.8)+
      scale_color_manual(values = c("sig"="red","nosig"="blue"),
                         breaks=c("sig","nosig"),name ="Significance",labels=c("p<0.05", "p>0.05"))+
      scale_size_area(max_size = 8)+
      geom_text_repel(label=aa$r)+labs(x=NULL,y="r")+
      geom_hline(yintercept=0,colour="grey50",linetype="dashed")+
      theme_classic()
  }
  return(p)
}


#' Performs graph-based permutation tests
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param metadata metadata
#' @param group one group name in columns of metadata
#' @param nperm numbers of permutations to perform
#' @param ... additional
#'
#' @return ggplots
#' @export
#'
#' @examples
#' data(otutab)
#' gp_dis_density(otutab,metadata,"Group")
grap_p_test<-function(otutab,metadata,group='Group', nperm = 999,...){
  lib_ps("igraph","phyloseq","phyloseqGraphTest","ggnetwork","intergraph")

  otumat = otutab
  #random taxon matrix for now. will update to real one later
  taxmat = matrix(sample(letters, 7*nrow(otumat), replace = TRUE), nrow = nrow(otumat), ncol = 7)
  rownames(taxmat) <- rownames(otumat)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  sampledata = sample_data(metadata)
  ps.expo = phyloseq(otu_table(otutab, taxa_are_rows =TRUE), tax_table(taxmat), sampledata)

  gt1 = graph_perm_test(ps.expo, group,nperm=nperm,...)
  cat('p-value:',gt1$pval,'\n')
  p1<-plot_test_network(gt1)
  p2<-plot_permutations(gt1)+
    labs(subtitle = paste0('p-value = ',gt1$pval,',   ',nperm,' permutations'))

  detach('package:phyloseqGraphTest')
  detach('package:phyloseq')
  detach('package:ggnetwork')
  detach('package:intergraph')
  return(list(p1,p2))
}


#' Group inter-intra density
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param group group vector
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(otutab)
#' gp_dis_density(otutab,metadata$Group)
gp_dis_density<-function(otutab,group){
  vegdist(t(otutab),method = "bray")%>%
    as.b_dist(group_df = group)%>%plot(.,mode=2)
}


#' RDA
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param env environmental factors
#' @param choose_var should choose variables? use forward step
#' @param norm should normalize?(default:T)
#' @param scale should scale species?(default:F)
#'
#' @return rda/cca
#' @export
#' @seealso \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @examples
#' data(otutab)
#' env=metadata[,12:16]
#' #RDA
#' myRDA(otutab,env)->phy.rda
#' RDA_plot(phy.rda,metadata$Group)
myRDA<-function(otutab,env,choose_var=F,norm=T,scale=F){
  lib_ps("vegan")
  data.frame(t(otutab))->dat
  if(norm)dat.h <- decostand(dat, "hellinger")else dat.h<-dat

  print(vegan::decorana(dat.h)->dca)
  cat('DCA分析，根据Axis lengths行的第一个值选择排序分析模型
  Axis Lengths >4.0-CCA(基于单峰模型，典范对应分析)；
  如果在3.0-4.0之间-RDA/CCA均可；
  如果小于3.0-RDA(基于线性模型，冗余分析)\n')
  (dcap<-dca$rproj%>%data.frame()%>%
      ggplot(.,aes(DCA1,DCA2))+labs(title = "DCA")+
      geom_point(size = 2)+ #可在这里修改点的透明度、大小
      geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
      geom_hline(yintercept = 0, color = 'gray', size = 0.4)+theme_classic())

  #env<-decostand(env,method = 'log',MARGIN = 2)#把环境因子进行log转化，以减少同一种环境因子之间本身数值大小造成的影响。

  #RDA
  (phy.rda=rda(dat.h~.,env,scale = scale))
  print('===============初始模型================')
  #print(anova(phy.rda, permutations = how(nperm = 999),by = 'terms'))
  #adonis2(dat.h~.,env,permutations = 9999, method="bray",by = 'terms')
  print('初始rda,vif>20表示共线性严重:')#vif>20表示共线性严重。
  print(vif.cca(phy.rda))
  cat('初始模型R2:',(R2a.all <- RsquareAdj(phy.rda)$adj.r.squared),'\n')

  #变量的选择
  if (choose_var) {
    # Compare the variance inflation factors
    spe.rda.all <- rda(dat.h ~ ., data = env)
    # Forward selection using forward.sel()
    #forward.sel(dat.h, env, adjR2thresh = R2a.all)#不能有非数值变量
    #或者使用ordistep
    mod0 <- rda(dat.h ~ 1, data = env)
    step.forward <-
      ordistep(mod0,
               scope = formula(spe.rda.all),
               direction = "forward",
               permutations = how(nperm = 499)
      )
    ## Parsimonious RDA，简化后的RDA
    print('==================筛选模型==============')
    #anova(step.forward, permutations = how(nperm = 999),by ='terms' )
    #adonis2(step.forward$call$formula,data = env,permutations = 999, method="bray")
    cat('筛选rda,vif>20表示共线性严重:\n')#vif>20表示共线性严重。
    print(vif.cca(step.forward))
    cat('筛选模型R2:',(R2a.pars <- RsquareAdj(step.forward)$adj.r.squared),'\n')
    phy.rda=step.forward
  }

  print('=============统计===========')
  B.sum=summary(phy.rda)
  cat(B.sum$constr.chi/B.sum$tot.chi,'constrained表示环境因子对群落结构差异的解释度','\n')
  cat(B.sum$unconst.chi/B.sum$tot.chi,'unconstrained表示环境因子对群落结构不能解释的部分\n')
  return(phy.rda)
}
#' @export
#' @rdname myRDA
myCCA<-function(otutab,env,choose_var=F,norm=T,scale=F){
  lib_ps("vegan")
  data.frame(t(otutab))->dat
  if(norm)dat.h <- decostand(dat, "hellinger")else dat.h<-dat

  print(decorana(dat.h)->dca)
  cat('DCA分析，根据Axis lengths行的第一个值选择排序分析模型
  Axis Lengths >4.0-CCA(基于单峰模型，典范对应分析)；
  如果在3.0-4.0之间-RDA/CCA均可；
  如果小于3.0-RDA(基于线性模型，冗余分析)\n')
  (dcap<-dca$rproj%>%data.frame()%>%ggplot(.,aes(DCA1,DCA2))+
      geom_point(size = 2)+ #可在这里修改点的透明度、大小
      geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
      geom_hline(yintercept = 0, color = 'gray', size = 0.4)+theme_classic())

  #env<-decostand(env,method = 'log',MARGIN = 2)#把环境因子进行log转化，以减少同一种环境因子之间本身数值大小造成的影响。

  #CCA
  (phy.cca=cca(dat.h~.,env))
  print('===============初始模型================')
  #print(anova(phy.cca, permutations = how(nperm = 999),by = 'terms'))
  #adonis2(dat.h~.,env,permutations = 9999, method="bray",by = 'terms')
  print('初始cca,vif>20表示共线性严重:')#vif>20表示共线性严重。
  print(vif.cca(phy.cca))
  cat('初始模型R2:',(R2a.all <- RsquareAdj(phy.cca)$adj.r.squared),'\n')
  #变量的选择
  if (choose_var) {
    # Compare the variance inflation factors
    spe.cca.all <- cca(dat.h ~ ., data = env)

    #print(R2a.all <- RsquareAdj(spe.cca.all)$adj.r.squared)
    # Forward selection using forward.sel()
    #forward.sel(dat.h, env, adjR2thresh = R2a.all)#不能有非数值变量
    #或者使用ordistep
    mod0 <- cca(dat.h ~ 1, data = env)
    step.forward <-
      ordistep(mod0,
               scope = formula(spe.cca.all),
               direction = "forward",
               permutations = how(nperm = 499)
      )
    ## Parsimonious cca，简化后的cca
    print('==================筛选模型==============')
    #anova(step.forward, permutations = how(nperm = 999),by ='terms' )
    #adonis2(step.forward$call$formula,data = env,permutations = 999, method="bray")
    cat('筛选cca,vif>20表示共线性严重:\n')#vif>20表示共线性严重。
    print(vif.cca(step.forward))
    cat('筛选模型R2:',(R2a.pars <- RsquareAdj(step.forward)$adj.r.squared),'\n')
    phy.cca=step.forward
  }
  print('=============统计===========')
  B.sum=summary(phy.cca)
  cat(B.sum$constr.chi/B.sum$tot.chi,'constrained表示环境因子对群落结构差异的解释度','\n')
  cat(B.sum$unconst.chi/B.sum$tot.chi,'unconstrained表示环境因子对群落结构不能解释的部分\n')
  return(phy.cca)
}
#' @param dist dist method for CAP, see \code{\link[vegan]{vegdist}}
#' @export
#' @rdname myRDA
myCAP<-function(otutab,env,choose_var=F,norm=T,scale=F,dist="bray"){
  lib_ps("vegan")
  data.frame(t(otutab))->dat
  if(norm)dat.h <- decostand(dat, "hellinger")else dat.h<-dat
  (phy.cap=capscale(dat.h~.,env,scale = scale,dist=dist ))
  print('===============初始模型================')
  #print(anova(phy.cap, permutations = how(nperm = 999),by = 'terms'))
  #adonis2(dat.h~.,env,permutations = 9999, method="bray",by = 'terms')
  print('初始cca,vif>20表示共线性严重:')#vif>20表示共线性严重。
  print(vif.cca(phy.cap))
  cat('初始模型R2:',(R2a.all <- RsquareAdj(phy.cap)$adj.r.squared),'\n')
  return(phy.cap)
}

#' Plot RDA res
#'
#' @param phy.rda rda/cca object
#' @param Group group vector to color
#' @param mode which mode choosen to plot
#' @param tri plot species arrows?
#' @param Topn how many variables to show?
#' @param rate segments length rate
#' @param sample_label show the samples label?
#' @param pal colors for group
#'
#' @seealso  \code{\link{myRDA}}
#' @return ggplot
#' @export
RDA_plot<-function(phy.rda,Group,mode=1,tri=F,Topn=10,rate=1,pal=NULL,sample_label=T){
  lib_ps("ggplot2","dplyr")
  if(is.null(pal)&!is.numeric(Group))pal=pcutils::get_cols(n = length(unique(Group)),pal = RColorBrewer::brewer.pal(5,"Set2"))
  getplotdat<-\(phy.rda,scale=1){
    #提取样方和环境因子排序坐标，前两轴，I 型标尺
    rda.scaling1 <- summary(phy.rda, scaling = scale)
    rda.site <- data.frame(rda.scaling1$sites)[1:2]
    rda.site$sample <- rownames(rda.site)
    rda.env <- data.frame(rda.scaling1$biplot)[1:2]
    rda.env$sample <- rownames(rda.env)
    rda.spe<-data.frame(rda.scaling1$species)[1:2]
    rda.spe$sample <- rownames(rda.spe)
    return(list(rda.site,rda.env,rda.spe))
  }

  plotdat=getplotdat(phy.rda)

  rda_eig <- (phy.rda$CCA$eig/sum(phy.rda$CCA$eig))[1:2]
  names(rda_eig)->approach

  colnames(plotdat[[2]])[1:2]<-c("RDA1","RDA2")
  colnames(plotdat[[3]])[1:2]<-c("RDA1","RDA2")

  plotdat1<-data.frame(plotdat[[1]],level=Group)
  colnames(plotdat1)[1:2]<-c("x1","x2")
  p<-plot_b_like(plotdat1,mode = mode,pal = pal,sample_label = sample_label)
  #envs
  p<-p+labs(x = paste(approach[1],': ', round(100 * rda_eig[1], 2), '%'),
         y = paste(approach[2],': ', round(100 * rda_eig[2], 2), '%'))+
    geom_segment(data = plotdat[[2]], aes(x = 0, y = 0, xend = RDA1,yend = RDA2),
                 arrow = arrow(length = unit(0.2, 'cm')), size = 0.5, color = 'blue') +
    geom_text(data = plotdat[[2]], aes(RDA1 * 1.1, RDA2 * 1.1, label = sample),
              color = 'blue', size = 4)

  if (!tri)return(p)

  if (tri){
    lib_ps("ggnewscale")
    plotdat[[3]]%>%
      transmute(contri=((RDA1^2/rda_eig[1])+(RDA2^2/rda_eig[2]))*100/2)%>%
      cbind(plotdat[[3]],.)->var_contri
    var_contri[,1:2]*rate->var_contri[,1:2]
    var_contri%>%top_n(Topn,contri)->var_contri
      p <- p+ggnewscale::new_scale_color()+
        geom_segment(data = var_contri, aes(x = 0, y = 0, xend = RDA1,yend = RDA2,color=contri),
                     size = 0.3,arrow = arrow(length=unit(0.1, "inches"))) +
        geom_text(data = var_contri, aes(x =  RDA1*1.05, y = RDA2*1.05,
                                      color=contri, label = sample), size = 2)+
        scale_color_gradientn(name="contribution",colours = c("#00AFBB", "#E7B800", "#FC4E07"))
  }
  return(p)
}


#' Envfit test for RDA result
#'
#' @param phy.rda a rda result
#' @param env environmental factors
#'
#' @return g_test object
#' @export
#' @seealso \code{\link[vegan]{envfit}}
#' @examples
#' data(otutab)
#' env=metadata[,12:16]
#' #RDA
#' myRDA(otutab,env)->phy.rda
#' enviftt(phy.rda,env)
#'
envfitt<-function(phy.rda,env,...){
  B.ef=envfit(phy.rda,env,...)#是做每一个环境因子与群落结构差异的相关性(解释量)
  cor_com<-data.frame(group=names(B.ef$vectors$r),
                      r2=B.ef$vectors$r,p_value=B.ef$vectors$pvals)
  cor_com$sig=cor_com$p_value<0.05
  class(cor_com)<-c("g_test","data.frame")
  return(cor_com)
}


#' Multi-table test with env
#'
#' @param g_otutab multi-otutabs with first column is group
#' @param env environmental factors
#'
#' @return a mant_g object
#' @export
#'
#' @examples
#' data(otutab)
#'cbind(group=rep(c('a','b','c'),c(200,100,192)),otutab)->g_otutab
#'metadata[,12:16,drop=F]->env
#'m_group_env(g_otutab,env)->mant_g
#'plot(mant_g,env)
m_group_env<-function(g_otutab,env){
  groups<-g_otutab$group%>%unique()
  all<-data.frame()
  for (i in groups){
    filter(g_otutab,group==i)[,-1]->tmp
    suppressWarnings(permanova(tmp,env,method = "mantel")->res)
    all=rbind(all,data.frame(spec=i,res))
  }
  all%>%mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                            labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                   pd = cut(p_value, breaks = c(-Inf, 0.01, 0.05, Inf),
                            labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))->all
  class(all)<-c("mant_g",class(all))
  return(all)
}


#' Plot mant_g object
#'
#' @return a ggplot
#' @exportS3Method
#'
#' @rdname m_group_env
plot.mant_g<-function(mant_g,env){
  lib_ps("ggcor")
  set_scale(c("#6D9EC1", "white", "#E46726"),type = "gradient2n")
  corp<-quickcor(env, type = "upper") +
    #geom_square() +
    geom_pie2()+
    anno_link(aes(colour = pd, size = rd), data = mant_g) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A288")) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))
  return(corp)
}


