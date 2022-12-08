
#' Calulate distance for otutab
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...
#' @param spe_nwk a phylo tree
#'
#' @return dist
#' @export
#' @seealso \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @examples
#' data(otutab)
#' b_dist(otutab)
#' makeNewick(taxonomy)->spe_nwk
#' b_dist(otutab)
b_dist<-function(otutab,method="bray",spe_nwk=NULL){
  totu=t(otutab)
  if(is.null(spe_nwk)){lib_ps("vegan");vegan::vegdist(totu,method = method)->a}
  else {
    lib_ps("picante")
    match.phylo.comm(spe_nwk,totu)->match_p
    match_p$phy->spe_nwk;match_p$comm->totu
    a<-switch (method,
               "unifrac" = picante::unifrac(totu, spe_nwk),
               "b_mpd" = picante::comdist(totu,cophenetic(spe_nwk)),
               "b_mntd"=picante::comdistnt(totu,cophenetic(spe_nwk)),
               "phylosor"=picante::phylosor(totu, spe_nwk)
               #"pcd"= picante::pcd(totu, spe_nwk),
               #picante::psd(samp,test)
               #picante::raoD(samp,test)
    )
  }
  return(a)
}

#' Transfer dist to b_dist
#'
#' @param dist a dist object
#' @param group_df a dataframe with rowname same to dist and one group column
#'
#' @return a b_dist with annotation by group
#' @export
#' @examples
#' data(otutab)
#' b_dist(otutab)%>%as.b_dist(.,group_df = metadata[,"Group",drop=F])->aa
#' plot(aa)
#' plot(aa,mode=2)
as.b_dist<-function(dist,group_df=NULL){
  #将dist矩阵转化为group注释的b_dist对象
  stopifnot(class(dist)=="dist")
  NST::dist.3col(dist)->aa
  if(!is.null(group_df)){
    group=group_df%>%unlist()
    if(F){
      #另一种方法
      as.matrix(a)->b
      group<-droplevels.factor(group)
      combn(levels(factor(group)),2)%>%apply(.,2,FUN = function(x)paste0(x[1],"_",x[2]))->pair
      list()->b_ls
      for (i in 1:nrow(b)){
        for (j in 1:nrow(b)){
          b_ls[[paste0(group[i],"_",group[j])]]<-c(b_ls[[paste0(group[i],"_",group[j])]],b[i,j])
        }
      }
      b_ls<-b_ls[pair]
      suppressMessages(as.data.frame(b_ls)%>%melt()->aa)
      colnames(aa)[2]<-"dis"
    }
    com<-function(group1,group2,levels){
      factor(c(group1,group2),levels = levels)%>%sort%>%paste(.,collapse = "_")
    }
    if (nrow(as.matrix(dist))!=length(group)) stop('group length wrong')
    aa<-left_join(aa,data.frame(name1=rownames(group_df),group1=group))
    aa<-left_join(aa,data.frame(name2=rownames(group_df),group2=group))
    aa%>%mutate(a=com(group1,group2,levels(group)))

    aa$variable<-apply(aa,1,function(x)com(x[4],x[5],levels(group)))
    aa%>%mutate(group=ifelse(group1==group2,"intra","inter"))->aa
    aa%>%mutate(variable=ifelse(group=="inter",variable,as.character(group1)))->aa
  }
  class(aa)<-c("b_dist",class(aa))
  return(aa)
}

#' Transfer b_dist to dist
#'
#' @param aa a b_dist object
#' @import reshape2
#' @return dist
#' @exportS3Method
as.dist.b_dist<-function(aa){
  lib_ps("reshape2","tibble")
  dcast(aa,name1~name2,value.var = "dis")%>%column_to_rownames("name1")->tmp
  rbind(rep(NA,ncol(tmp)),tmp)->tmp
  rownames(tmp)[1]= setdiff(colnames(tmp),rownames(tmp))
  cbind(tmp,rep(NA,nrow(tmp)))->tmp
  colnames(tmp)[ncol(tmp)]=setdiff(rownames(tmp),colnames(tmp))
  tmp[rownames(tmp),rownames(tmp)]->tmp
  tmp[is.na(tmp)]<-0
  tmp+t(tmp)->tmp
  return(as.dist(tmp))
}

#' Plot dist
#'
#' @param aa a dist
#' @param group_df a dataframe with rowname same to dist and one group column
#' @param ... additional
#'
#' @return a pheatmap
#' @exportS3Method
#' @rdname as.b_dist
plot.dist<-function(aa,group_df=NULL,...){
  lib_ps("pheatmap","dplyr")
  ab<-as.matrix(aa)%>%as.data.frame()
  if(!is.null(group_df)){
    group=group_df%>%unlist()
    ann=group_df%>%arrange(factor(group))
    ab=ab[rownames(ann),rownames(ann)]
    pheatmap(ab,annotation_row = ann,annotation_col = ann,cluster_rows = F,
             cluster_cols = F,border_color = F,...)
  }
  else{
    pheatmap(ab)
  }
}

#' Plot b_dist
#'
#' @param aa a b_dist
#' @param mode 1~3
#' @param c_group "inter" or "intra" or both to plot
#' @param ... additional
#'
#' @return a ggplot or pheatmap
#' @exportS3Method
#' @rdname as.b_dist
#'
plot.b_dist<-function(aa,mode=1,c_group="inter",...){
  lib_ps("dplyr")
  if(mode==1){
    aa$variable=factor(aa$variable,levels = unique(c(levels(aa$group1),unique(aa$variable))))
    filter(aa,group%in%c_group)->aa1
    p=aa1%>%group_box(.,3,group = .$variable,...)
    return(p)
  }
  if(mode==2){
    pl=ggplot(aa,aes(x=dis,col=group,fill=group))+
      geom_density(alpha=0.5)+
      geom_rug()+
      geom_vline(aes(xintercept=median,col=group),data = aa%>%group_by(group)%>%summarise(median=median(dis)),linetype=2)+
      scale_color_manual(values = c( "#E88493","#51A4E0"))+
      scale_fill_manual(values = c( "#E88493","#51A4E0"))+ labs(x="Bray-Curtis Distance")
    #获得ggplot绝对画板大小
    (wilcox.test(aa$dis~aa$group)->wt)
    ggplot_build(pl)->plot
    p=pl+annotate('text',x = 0.9*max(plot$data[[1]]$x),y = max(plot$data[[1]]$y),label=paste0('p=',signif(wt$p.value,4)))
    return(p)
  }
  if(mode==3){
    ann=data.frame(name=c(aa$name1,aa$name2),group=c(aa$group1,aa$group2))%>%distinct()%>%
      column_to_rownames("name")
    plot(as.dist(aa),ann)
  }
}

#Stochastic or Determine=========
#' Calculate NST for each group
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param group_df a dataframe with rowname and one group column
#' @param threads default:4
#' @param file filename to save
#' @param rep repeat numbers: suggest 999
#' @return a b_dist object, dis is MSTij
#' @export
#'
#' @references 1. Ning, D., Deng, Y., Tiedje, J. M. & Zhou, J. A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892–16898 (2019).
#' @examples
#' data(otutab)
#' nst(otutab,metadata[,"Group",drop=F])->nst_res
#' plot(nst_res,c_group = "intra")+geom_hline(yintercept = 0.5,lty=2)
nst<-function(otutab,group_df,threads=4,file=NULL,rep=20){
  lib_ps("NST","dplyr")
  tnst=NST::tNST(comm=t(otutab), group=group_df,rand=rep,output.rand=TRUE, nworker=threads)
  if(is.null(file))file="tsnt_res"
  tnst$index.pair%>%mutate(dis=MST.ij.ruzicka)%>%as.dist.b_dist()->NST_ij
  as.b_dist(NST_ij,group_df)->NST
  save(tnst,NST,file = paste0(file,".rda"))
  return(NST)
}

#' Sloan Neutral Model
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param model fit method, one of "nls","mle"
#'
#' @return ncm_res
#' @export
#' @references 1. Sloan, W. T. et al. Quantifying the roles of immigration and chance in shaping prokaryote community structure. Environmental Microbiology 8, 732–740 (2006).
#' @examples
#' data(otutab)
#' ncm(otutab)->ncm_res
#' plot(ncm_res)
ncm<-function(otutab,model="nls"){
  lib_ps("Hmisc")
  otutab<-otutab[rowSums(otutab)>0,]
  spp<-t(otutab)
  otu.pa <- (spp > 0) * 1
  freq <- apply(otu.pa, 2, mean)
  N <- mean(apply(spp, 1, sum))
  p <- apply(spp, 2, function(x) mean(x))/N
  d = 1/N
  if(model=="nls"){
    lib_ps("minpack.lm")
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- minpack.lm::nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))

    #coef(m.fit) #get the m value
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
    pred.ci <- Hmisc::binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])

    stats=data.frame(m=coef(m.fit),R2=Rsqr,N=N,samples=ncol(otutab),tax=nrow(otutab),d=d)}

  else if(model=="mle"){
    lib_ps("bbmle")
    #Maximum Likelihood Estimation
    neutral.ll <- function(m, sigma) {
      R = freq - pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
      -sum(dnorm(R, 0, sigma, log = TRUE))
    }
    m.mle <- bbmle::mle2(neutral.ll, start = list(m = 0.01, sigma = 0.1),
                         method = "Nelder-Mead")

    freq.pred <- pbeta(d, N * m.mle@coef["m"] * p, N * m.mle@coef["m"] * (1 - p), lower.tail = FALSE)
    pred.ci <- Hmisc::binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

    gRsqr <- 1 - exp(-as.numeric(logLik(m.mle))/length(p))
    bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])
    stats=data.frame(m=m.mle@coef["m"], R2= gRsqr,N=N,samples=ncol(otutab),tax=nrow(otutab),d=d)
    }
  bacnlsALL$group = NA
  bacnlsALL$group[bacnlsALL[,2]<bacnlsALL[,4]]="Below" ##低于下界
  bacnlsALL$group[bacnlsALL[,2]>bacnlsALL[,5]]="Above" ##高于上界
  bacnlsALL$group[(bacnlsALL[,2]>=bacnlsALL[,4])&(bacnlsALL[,2]<=bacnlsALL[,5])]="In"##中间
  bacnlsALL<-bacnlsALL[rownames(otutab),]
  res=list(stats=stats,ncm_data=bacnlsALL)
  class(res)<-"ncm_res"
  return(res)
}

#' Plot ncm_res
#'
#' @param res a ncm_res object
#'
#' @return ggplot
#' @exportS3Method
#'
#' @seealso \code{\link{ncm}}
plot.ncm_res<-function(ncm_res){
  lib_ps("ggpubr","patchwork")
  out = ncm_res[[2]]
  lincol='#4a8fb8'
  p1 = ggplot() +
    geom_line(data = out,aes(x=log(p),y=freq.pred),linewidth = 1.2,linetype = 1,col=lincol)+
    geom_line(data = out,aes(x=log(p),y=Lower),linewidth = 1.2,linetype = 2,col=lincol)+
    geom_line(data = out,aes(x=log(p),y=Upper),linewidth = 1.2,linetype = 2,col=lincol)+
    xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")

  mycols<-c("Above"="#069870","Below"="#e29e02","In"="#1e353a")

  p2 = p1 + geom_point(data = out,aes(x=log(p),y=freq,color = group),size = 1)+
    scale_colour_manual(values = mycols)+
    annotate("text",x=log(max(out$p))-2,y=0.05,label=paste("Nm = ",sprintf("%.0f",ncm_res[[1]][1]*ncm_res[[1]][3]),sep=''),size=4)+
    annotate("text",x=log(max(out$p))-2,y=0.1,label=paste0('R2 = ',round(ncm_res[[1]][2],3)),size=4)+
    theme_pubr()+theme(legend.position = c(0.85,0.3),
                       legend.title = element_blank(),legend.background=element_rect(I(0)))

  out$group%>%table()%>%as.data.frame()->ad
  colnames(ad)<-c("type","n")
  ad$fraction = ad$n / sum(ad$n)
  ad$ymax = cumsum(ad$fraction)
  ad$ymin = c(0, head(ad$ymax, n = -1))
  ad$rate_per<-paste(as.character(round(100*ad$fraction,1)),'%',sep='')
  pie=ggplot(data = ad, aes(fill = type, ymax = ymax, ymin = ymin, xmax = 3, xmin = 1.5)) +
    geom_rect(show.legend = F,alpha=0.8) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", title = "",fill='') +
    xlim(c(0, 3)) +theme_light() +
    theme(panel.grid=element_blank()) + ## 去掉白色外框
    theme(axis.text=element_blank()) + ## 把图旁边的标签去掉
    theme(axis.ticks=element_blank()) + ## 去掉左上角的坐标刻度线
    theme(panel.border=element_blank()) + ## 去掉最外层的正方形边框
    #geom_text(aes(x = 3.6, y = ((ymin+ymax)/2),label = type) ,size=2)+
    geom_text(aes(x = 2.2, y = ((ymin+ymax)/2),label = rate_per) ,size=2,col="white")+
    scale_fill_manual(values = mycols)+
    theme(plot.background=element_rect(I(0),linetype=0),
          panel.background=element_rect(I(0)))#去除图片背景白色

  p2+ inset_element(pie,left = -0.1,bottom = 0.6,right = 0.3,top = 1.1)
}


#' Calculate beta_NTI
#'
#' @param phylo a phylo object
#' @param otutab otutab
#' @param beta.reps how many simulation performed?
#' @param threads use how many threads to calculate (default:4)
#'
#' @return a dist: b_NTI
#' @export
#'
#' @examples
#' data(otutab)
#' makeNewick(taxonomy)->phylo
#' b_NTI(phylo,otutab)->bnti_res
#' data.frame(type=c("Stochastic","Homo_S","Heter_S"),
#'           number=c(sum(abs(bnti_res)<2),
#'              sum(bnti_res<(-2)),sum(bnti_res>2)))->com_pro
#' pcutils::gghuan(com_pro)
#'
b_NTI<-function(phylo,otutab,beta.reps=9,weighted=T,threads=4){
  lib_ps("picante")
  #match tree and otutab (important)
  phy_otu_m = match.phylo.data(phylo, otutab)
  #comdistnt:Calculates inter-community mean nearest taxon distance
  #cophenetic:Cophenetic Distances for a Hierarchical Clustering

  beta.mntd.weighted = as.matrix(comdistnt(t(phy_otu_m$data),
                                           cophenetic(phy_otu_m$phy),abundance.weighted=weighted))

  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(phy_otu_m$data), ncol(phy_otu_m$data),beta.reps))
if(threads==1){
  #no parallel
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(
      comdistnt(t(phy_otu_m$data),taxaShuffle(cophenetic(phy_otu_m$phy)),
                abundance.weighted=weighted,exclude.conspecifics = F))
    #print(c(date(),rep))
  }
}
else if(threads>1){
  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =beta.reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  rand.weighted.bMNTD.comp <- foreach (rep = 1:beta.reps,
                                       .options.snow = opts,
                                       .packages = c("picante")) %dopar%{
    as.matrix(comdistnt(t(phy_otu_m$data),
                        taxaShuffle(cophenetic(phy_otu_m$phy)),
                        abundance.weighted=T,exclude.conspecifics = F))
  }
  stopCluster(cl)
  rand.weighted.bMNTD.comp<-simplify2array(rand.weighted.bMNTD.comp)
}
  weighted.bNTI = matrix(c(NA),nrow=ncol(phy_otu_m$data),ncol=ncol(phy_otu_m$data))
  for(columns in 1:(ncol(phy_otu_m$data)-1)) {
    for(rows in (columns+1):ncol(phy_otu_m$data)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] -mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")
    }
  }
  rownames(weighted.bNTI) = colnames(phy_otu_m$data);
  colnames(weighted.bNTI) = colnames(phy_otu_m$data);
  return(as.dist(weighted.bNTI))
}


#' Calculate RCbray-curtis
#'
#' @param otutab otutab
#' @param reps how many simulation performed?
#' @param threads use how many threads to calculate (default:4)
#' @param classic_metric standardizes the metric to range from -1 to 1
#' @param split_ties adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended
#' @details Parallelized version of the Raup-Crick algorithm for "abundance" data (Stegen et al. 2013).
#' @return a dist
#' @export
#'
#' @examples
#' data(otutab)
#' makeNewick(taxonomy)->phylo
#' b_NTI(phylo,otutab)->bnti_res
#' RCbray(otutab,reps=9)->rc_res
#'
#' data.frame(type=factor(c("Homo_S","Heter_S","Homo_D","D_limit","Undominated"),
#'              levels = c("Homo_S","Heter_S","Homo_D","D_limit","Undominated")),
#'            number=c(sum(bnti_res<(-2)),sum(bnti_res>2),
#'              sum((abs(bnti_res)<2)&(abs(rc_res)<0.95)),
#'              sum((abs(bnti_res)<2)&(rc_res<(-0.95))),
#'              sum((abs(bnti_res)<2)&(rc_res>0.95))))->com_pro
#' gghuan(com_pro,reorder = F)
RCbray <- function(otutab, reps=9, threads=4, classic_metric=T, split_ties=TRUE){
  lib_ps("vegan")
  lib_ps("parallel")
  lib_ps("doSNOW")
  com=t(otutab)
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  bray.rand <- foreach(randomize = 1:reps,.options.snow = opts,.packages = c("vegan")) %dopar% {
               null.dist <- com*0
               for(i in 1:nrow(com)){
                 com.pa <- (com>0)*1
                 gamma<-ncol(com)
                 occur<-apply(com>0, MARGIN=2, FUN=sum)
                 abundance<-apply(com, MARGIN=2, FUN=sum)
                 com1 <- rep(0,gamma)
                 com1[sample(1:gamma, sum(com.pa[i,]), replace=FALSE, prob=occur)]<-1
                 com1.samp.sp = sample(which(com1>0), (sum(com[i,])-sum(com1)),
                                       replace=TRUE,prob=abundance[which(com1>0)]);
                 com1.samp.sp = cbind(com1.samp.sp,1)
                 com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum))
                 colnames(com1.sp.counts) = 'counts'
                 com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
                 com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
                 x <- com1
                 null.dist[i,] <- x
                 rm('com1.samp.sp','com1.sp.counts')
               }
               as.matrix(vegdist(null.dist, "bray"))
  }

  stopCluster(cl)
  ## Calculate beta-diversity for obs metacommunity
  bray.obs <- as.matrix(vegdist(com, "bray"))

  ##how many null observations is the observed value tied with?
  null_bray_curtis <- bray.rand
  num_exact_matching_in_null <- lapply(null_bray_curtis, function(x) x==bray.obs)
  num_exact_matching_in_null <- apply(simplify2array(num_exact_matching_in_null), 1:2, sum)

  ##how many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null <- lapply(null_bray_curtis, function(x) (x<bray.obs)*1)
  num_less_than_in_null <- apply(simplify2array(num_less_than_in_null), 1:2, sum)

  rc = (num_less_than_in_null)/reps; # rc;

  if(split_ties){
    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
  };
  if(!classic_metric){
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    rc = (rc-.5)*2
  };
  return(as.dist(rc))
}


#=========iCAMP====
#里面有更多的写好的函数
if(F){
  #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247487267&idx=1&sn=6b776e219d1a782adc96b3f4937c2d7a&chksm=ea71f1e8dd0678fe0f5e3e03a06f98483cf2923c0df0f600c821dcda8c109217924afbe91ad4&scene=21#wechat_redirect
  #wd0=getwd() # please change to the folder you want to save the pd.big output.
  pd.wd=paste0(getwd(),"/pdbig.icampbig")
  nworker=4 # parallel computing thread number
  rand.time=20 # usually use 1000 for real data.
  bin.size.limit=5 # for real data, usually use a proper number
  library(iCAMP)
  icamp.out=icamp.big(comm=t(otutab),tree=phylo,pd.wd=pd.wd,
                      rand=rand.time, nworker=nworker,
                      bin.size.limit=bin.size.limit)


  #计算树的距离矩阵
  pd=iCAMP::pdist.p(phylo,nworker = 4)
  pd.big=iCAMP::pdist.big(phylo)
  #计算bNTI
  iCAMP::bmpd(t(otutab),pd)
  bnti=iCAMP::bNTIn.p(t(otutab),pd,rand = 20)
  #计算RC
  iCAMP::RC.pc(t(otutab),rand = 20)
  #计算生态位差异
  iCAMP::dniche()
}

#beta-diversity decomposition
if(F){
#https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247484144&idx=1&sn=a181f9d6a5e47681727901b80907c24e&chksm=ea71fc3bdd06752df3cac0865ac94736138abde1d85741625e4ccbc3c6d19e34896dce06b6cd&scene=21#wechat_redirect
lib_ps("betapart")
data(ceram.s)
#得到三个矩阵，两两成对比较。分别为总beta多样性jac or sor,物种转换jtu or sim,物种增减jne or sne.
ceram.core.s<-betapart.core(ceram.s)
ceram.dist.jac<-beta.pair(ceram.core.s, index.family="jac")
ceram.dist.sor<-beta.pair(ceram.core.s, index.family="sor")
ceram.multi.jac<-beta.multi(ceram.core.s, index.family="jac")
ceram.multi.sor<-beta.multi(ceram.core.s, index.family="sor")
}


#======ecological niche=====
##https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247485757&idx=1&sn=88678b4be64ac8774f0672b0155c295c&chksm=ea71f7f6dd067ee0407c9ee894a1d6ba8e052fa8a66e0df31e950309093f4549b0e3a5aa6830&scene=21#wechat_redirect
if(F){
  lib_ps("spaa")
  #### Niche width and niche overlap
  #niche.width计算生态位宽度
  #mat:列为物种，行为样本
  #method:计算方法
  niche_w<-niche.width(t(otutab), method = "levins")%>%t%>%as.data.frame()

  #每两个物种之间生态位的重叠
  niche.overlap(datasample[,1:3], method = "levins")
  #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247486549&idx=1&sn=05c283e4cf0afc4d338f1ca1c11dc5fc&chksm=ea71f29edd067b8829272ee565242ad9b08d18c4813d9e025144b41b9506c4cba40b754e2b6e&scene=21#wechat_redirect
  lib_ps("MicroNiche")
  ?MicroNiche::levins.Bn()
  ?MicroNiche::hurlberts.Bn()
  MicroNiche::hurlberts.Bn(otutab)
  ?MicroNiche::feinsingers.PS()
  #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247485087&idx=1&sn=95d121fc4080aa9af2a61de22694093d&chksm=ea71f854dd067142c0b201657c57aeb6102cef5f738fc89380783efa8c413580515b8536e5a6&scene=21#wechat_redirect
  lib_ps("indicspecies")

  #labdsv包的indval函数可计算群落中的指示物种
  lib_ps("labdsv")
  labdsv::indval()
  #分析微生物群落和环境因子相关性的工具vegan::bioenv
  vegan::bioenv()
  #表征微生物抵抗性的一个指标(要有实验-对照的关系才能算)
  #Rs=1-(2*|D0|/(C0+|D0|))
  #其中C0为对照组样本的α多样性指数，D0为对照组样本多样性和实验组样本α多样性之间的差值。
  #Rs取值为-1到1。越接近1表明处理效应非常小，即微生物抵抗性很强；Rs越小表明处理效应越强，即微生物抵抗力越弱。


}
devtools::install_github("Corentin-Gibert-Paleontology/DNCImper")

#========Cohesion======
#https://mp.weixin.qq.com/s?__biz=MzAwODk1Njk5MA==&mid=2247485401&idx=1&sn=30f4c324e1013482b98aa5b5b4f5cdc1&scene=21#wechat_redirect
