#Network========

#' Correlation network, species-interaction network for omics
#'
#' @return No value
#' @export
cor_net=function(){
  if(!requireNamespace("MetaNet")){
  message('All network related functions have been migrated to `MetaNet`, please install `MetaNet`:
remotes::install_github("Asa12138/MetaNet")
visit the website: https://github.com/Asa12138/MetaNet')
  }
  else {
    message('All network related functions have been migrated to `MetaNet`, you have already installed `MetaNet`, please:
library(MetaNet)
help(c_net_build)')
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
#'
#' @return a b_dist object, dis is MSTij
#' @export
#'
#' @references 1. Ning, D., Deng, Y., Tiedje, J. M. & Zhou, J. A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892–16898 (2019).
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' nst(otutab,metadata["Group"])->nst_res
#' plot(nst_res,c_group = "intra")+geom_hline(yintercept = 0.5,lty=2)+ylab("NST")
#' }}
nst<-function(otutab,group_df,threads=1,file=NULL,rep=20){
  lib_ps("NST","dplyr",library = FALSE)
  tnst=NST::tNST(comm=t(otutab), group=group_df,rand=rep,output.rand=TRUE, nworker=threads)

  tnst$index.pair%>%dplyr::mutate(dis=MST.ij.ruzicka)%>%as.dist.b_dist()->NST_ij

  as.b_dist(NST_ij,group_df)->NST

  if(is.null(file))file = paste0("nst_res_",date(),".RDS")
  saveRDS(tnst,file = file)
  message(paste0("Result saved as ",file))
  return(NST)
}

#' Calculate b_NTI and RC_bray for each group
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param phylo  a phylo object
#' @param group_df a dataframe with rowname and one group column
#' @param threads default:4
#' @param file filename to save
#' @param rep repeat numbers: suggest 999
#'
#' @return a b_dist object, dis is MSTij
#' @export
#'
#' @references 1. Ning, D., Deng, Y., Tiedje, J. M. & Zhou, J. A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892–16898 (2019).
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' df2tree(taxonomy)->phylo
#' nti_rc(otutab,phylo,metadata["Group"])->nti_res
#' plot(nti_res)
#' }}
nti_rc<-function(otutab,phylo,group_df,threads=1,file=NULL,rep=20){
  lib_ps("NST","dplyr",library = FALSE)
  pdist=stats::cophenetic(phylo)

  bnti_res<-NST::pNST(comm = t(otutab),pd = pdist, group = group_df,
                 pd.wd = tempdir(),rand = rep, nworker = threads, SES = TRUE, RC = TRUE)

  bnti_res$index.pair%>%dplyr::mutate(dis=bNTI.wt)%>%as.dist.b_dist()->b_NTI_ij
  as.b_dist(b_NTI_ij,group_df)->b_NTI
  bnti_res$index.pair%>%dplyr::mutate(dis=RC.bMNTD.wt)%>%as.dist.b_dist()->RC_ij
  as.b_dist(RC_ij,group_df)->RC

  b_NTI%>%dplyr::mutate(bNTI=dis,RC=RC$dis)%>%dplyr::filter(group=="intra")%>%dplyr::select(name1,name2,variable,bNTI,RC)->NTI_RC
  NTI_RC%>%dplyr::mutate(type=ifelse(bNTI<(-2),"Homo_S",
                              ifelse(bNTI>2,"Heter_S",
                                     ifelse(RC<(-0.95),"D_limit",
                                            ifelse(RC>0.95,"Undominated","Homo_D")))))->NTI_RC

  if(is.null(file))file = paste0("nti_rc_res_",date(),".RDS")
  saveRDS(bnti_res,file = file)
  message(paste0("Result saved as ",file))
  class(NTI_RC)=c("NTI_RC","data.frame")
  return(NTI_RC)
}

#' Plot NTI_RC object
#'
#' @param x NTI_RC object
#' @param ... pass to \code{\link[pcutils]{stackplot}}
#'
#' @return ggplot
#' @exportS3Method
#' @method plot NTI_RC
#' @rdname nti_rc
plot.NTI_RC=function(x,...){
  nti_res=x
  nti_res$type=factor(nti_res$type,levels = c("Homo_S","Heter_S","Homo_D","D_limit","Undominated"))
  table(nti_res$type,nti_res$variable)%>%reshape2::melt()%>%reshape2::acast(Var1~Var2)->com_p
  pcutils::stackplot(com_p,...)
}

#' Sloan Neutral Model
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param model fit method, one of "nls","mle"
#'
#' @return ncm_res
#' @export
#' @references 1. Sloan, W. TRUE. et al. Quantifying the roles of immigration and chance in shaping prokaryote community structure. Environmental Microbiology 8, 732–740 (2006).
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' ncm(otutab)->ncm_res
#' plot(ncm_res)
#' }}
ncm<-function(otutab,model="nls"){
  lib_ps("Hmisc",library = FALSE)
  otutab<-otutab[rowSums(otutab)>0,]
  spp<-t(otutab)
  otu.pa <- (spp > 0) * 1
  freq <- apply(otu.pa, 2, mean)
  N <- mean(apply(spp, 1, sum))
  p <- apply(spp, 2, function(x) mean(x))/N
  d = 1/N
  if(model=="nls"){
    lib_ps("minpack.lm",library = FALSE)
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- minpack.lm::nlsLM(freq ~ stats::pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))

    #coef(m.fit) #get the m value
    freq.pred <- stats::pbeta(d, N*stats::coef(m.fit)*p, N*stats::coef(m.fit)*(1 -p), lower.tail=FALSE)
    pred.ci <- Hmisc::binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])

    stats=data.frame(m=stats::coef(m.fit),R2=Rsqr,N=N,samples=ncol(otutab),tax=nrow(otutab),d=d)}

  else if(model=="mle"){
    lib_ps("bbmle",library = FALSE)
    #Maximum Likelihood Estimation
    neutral.ll <- function(m, sigma) {
      R = freq - stats::pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
      -sum(stats::dnorm(R, 0, sigma, log = TRUE))
    }
    m.mle <- bbmle::mle2(neutral.ll, start = list(m = 0.01, sigma = 0.1),
                         method = "Nelder-Mead")

    freq.pred <- stats::pbeta(d, N * m.mle@coef["m"] * p, N * m.mle@coef["m"] * (1 - p), lower.tail = FALSE)
    pred.ci <- Hmisc::binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

    gRsqr <- 1 - exp(-as.numeric(stats::logLik(m.mle))/length(p))
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
#' @param x a ncm_res object
#' @param ... add
#' @param mycols mycols
#' @param text_position text_position
#'
#' @return ggplot
#' @exportS3Method
#' @method plot ncm_res
#' @rdname ncm
plot.ncm_res<-function(x,mycols=c("Above"="#069870","Below"="#e29e02","In"="#1e353a"),text_position=NULL,...){
  ncm_res=x
  lib_ps("ggpubr","patchwork",library = FALSE)
  out = ncm_res[[2]]
  lincol='#4a8fb8'
  p1 = ggplot() +
    geom_line(data = out,aes(x=log(p),y=freq.pred),linewidth = 1.2,linetype = 1,col=lincol)+
    geom_line(data = out,aes(x=log(p),y=Lower),linewidth = 1.2,linetype = 2,col=lincol)+
    geom_line(data = out,aes(x=log(p),y=Upper),linewidth = 1.2,linetype = 2,col=lincol)+
    xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")

  if(is.null(text_position)){
    text_position=c(log(max(out$p))-2,0.05)
  }

  p2 = p1 +
    geom_point(data = out,aes(x=log(p),y=freq,color = group),size = 1)+
    scale_colour_manual(values = mycols)+
    annotate("text",x=text_position[1],y=text_position[2],label=paste("Nm = ",sprintf("%.0f",ncm_res[[1]][1]*ncm_res[[1]][3]),sep=''),size=4)+
    annotate("text",x=text_position[1],y=text_position[2]+0.1,label=paste0('R2 = ',round(ncm_res[[1]][2],3)),size=4)+
    pctax_theme+theme(legend.position = c(0.85,0.3),
                       legend.title = element_blank(),legend.background=element_rect(I(0)))

  out$group%>%table()%>%as.data.frame()->ad
  colnames(ad)<-c("type","n")
  pie=pcutils::gghuan(ad,name = FALSE,text_params = list(size=2.5))+
    xlim(0.2,3.3)+
    scale_fill_manual(values = mycols)+
    theme(plot.background=element_rect(I(0),linetype=0),
          panel.background=element_rect(I(0)))#去除图片背景白色

  p2+ patchwork::inset_element(pie,left = -0.1,bottom = 0.6,right = 0.3,top = 1.1)

}


#' Calculate beta_NTI
#'
#' @param phylo a phylo object
#' @param otutab otutab
#' @param beta.reps how many simulation performed?
#' @param threads use how many threads to calculate (default:4)
#' @param weighted logical
#'
#' @return a dist: b_NTI
#' @export
b_NTI1<-function(phylo,otutab,beta.reps=9,weighted=TRUE,threads=1){
  lib_ps("picante",library = FALSE)
  #match tree and otutab (important)
  phy_otu_m = picante::match.phylo.data(phylo, otutab)
  #picante::comdistnt:Calculates inter-community mean nearest taxon distance
  #cophenetic:Cophenetic Distances for a Hierarchical Clustering

  beta.mntd.weighted = as.matrix(picante::comdistnt(t(phy_otu_m$data),
                                           stats::cophenetic(phy_otu_m$phy),abundance.weighted=weighted))

  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(phy_otu_m$data), ncol(phy_otu_m$data),beta.reps))
if(threads==1){
  #no parallel
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(
      picante::comdistnt(t(phy_otu_m$data),picante::taxaShuffle(stats::cophenetic(phy_otu_m$phy)),
                abundance.weighted=weighted,exclude.conspecifics = FALSE))
    #print(c(date(),rep))
  }
}
else if(threads>1){
  #parallel
  pcutils::lib_ps("foreach","doSNOW","snow")
  pb <- utils::txtProgressBar(max =beta.reps, style = 3)
  opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
  cl <- snow::makeCluster(threads)
  doSNOW::registerDoSNOW(cl)
  rand.weighted.bMNTD.comp <- foreach::foreach (rep = 1:beta.reps,
                                       .options.snow = opts,
                                       .packages = c("picante")) %dopar%{
    as.matrix(picante::comdistnt(t(phy_otu_m$data),
                                 picante::taxaShuffle(stats::cophenetic(phy_otu_m$phy)),
                        abundance.weighted=TRUE,exclude.conspecifics = FALSE))
  }
  snow::stopCluster(cl)
  gc()
  pcutils::del_ps("doSNOW","snow","foreach")
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
#'
#' @details Parallelized version of the Raup-Crick algorithm for "abundance" data (Stegen et al. 2013).
#' @return a dist
#' @export
#'
#' @examples
#' \dontrun{
#' \donttest{
#' data(otutab,package = "pcutils")
#' df2tree(taxonomy)->phylo
#' b_NTI1(phylo,otutab)->bnti_res
#' RCbray1(otutab,reps=9)->rc_res
#'
#' data.frame(type=factor(c("Homo_S","Heter_S","Homo_D","D_limit","Undominated"),
#'              levels = c("Homo_S","Heter_S","Homo_D","D_limit","Undominated")),
#'            number=c(sum(bnti_res<(-2)),sum(bnti_res>2),
#'              sum((abs(bnti_res)<2)&(abs(rc_res)<0.95)),
#'              sum((abs(bnti_res)<2)&(rc_res<(-0.95))),
#'              sum((abs(bnti_res)<2)&(rc_res>0.95))))->com_pro
#' pcutils::gghuan(com_pro,reorder = FALSE)
#' }}
RCbray1 <- function(otutab, reps=9, threads=1, classic_metric=TRUE, split_ties=TRUE){
  lib_ps("vegan",library = FALSE)
  com=t(otutab)

  pcutils::lib_ps("foreach","doSNOW","snow")
  pb <- utils::txtProgressBar(max =reps, style = 3)
  opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
  cl <- snow::makeCluster(threads)
  doSNOW::registerDoSNOW(cl)

  bray.rand <- foreach::foreach(randomize = 1:reps,.options.snow = opts,.packages = c("vegan")) %dopar% {
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
               as.matrix(vegan::vegdist(null.dist, "bray"))
  }

  snow::stopCluster(cl)
  gc()
  pcutils::del_ps("doSNOW","snow","foreach")
  ## Calculate beta-diversity for obs metacommunity
  bray.obs <- as.matrix(vegan::vegdist(com, "bray"))

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
if(FALSE){
  #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247487267&idx=1&sn=6b776e219d1a782adc96b3f4937c2d7a&chksm=ea71f1e8dd0678fe0f5e3e03a06f98483cf2923c0df0f600c821dcda8c109217924afbe91ad4&scene=21#wechat_redirect
  #wd0=getwd() # please change to the folder you want to save the pd.big output.
  tax_tree<-df2tree(taxonomy)
  comm=t(otutab)
  tree=tax_tree
  picante::match.phylo.comm(tree,comm)
  # since need to save some output to a certain folder,
  # the following code is set as 'not test'.
  # but you may test the code on your computer
  # after change the folder path for 'pd.wd'.
  ## No test:
  wd0=getwd() # please change to the folder you want to save the pd.big output.
  pd.wd=paste0(tempdir(),"/pdbig.icampbig")
  nworker=4 # parallel computing thread number
  rand.time=99 # usually use 1000 for real data.

  bin.size.limit=48 # for real data, usually use a proper number
  # according to phylogenetic signal test or try some settings
  # then choose the reasonable stochasticity level.
  # our experience is 12, or 24, or 48.
  # but for this example dataset which is too small, have to use 5.

  icamp.out=iCAMP::icamp.big(comm=comm,tree=tree,pd.wd=pd.wd,
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
if(FALSE){
  # #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247484144&idx=1&sn=a181f9d6a5e47681727901b80907c24e&chksm=ea71fc3bdd06752df3cac0865ac94736138abde1d85741625e4ccbc3c6d19e34896dce06b6cd&scene=21#wechat_redirect
  # lib_ps("betapart",library = FALSE)
  # data(ceram.s)
  # #得到三个矩阵，两两成对比较。分别为总beta多样性jac or sor,物种转换jtu or sim,物种增减jne or sne.
  # ceram.core.s<-betapart::betapart.core(ceram.s)
  # ceram.dist.jac<-betapart::beta.pair(ceram.core.s, index.family="jac")
  # ceram.dist.sor<-betapart::beta.pair(ceram.core.s, index.family="sor")
  # ceram.multi.jac<-betapart::beta.multi(ceram.core.s, index.family="jac")
  # ceram.multi.sor<-betapart::beta.multi(ceram.core.s, index.family="sor")
}


#======ecological niche=====
##https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247485757&idx=1&sn=88678b4be64ac8774f0672b0155c295c&chksm=ea71f7f6dd067ee0407c9ee894a1d6ba8e052fa8a66e0df31e950309093f4549b0e3a5aa6830&scene=21#wechat_redirect

if(FALSE){
  # lib_ps("spaa",library = FALSE)
  # #### Niche width and niche overlap
  # #niche.width计算生态位宽度
  # #mat:列为物种，行为样本
  # #method:计算方法
  # niche_w<-spaa::niche.width(t(otutab), method = "levins")%>%t%>%as.data.frame()
  #
  # #每两个物种之间生态位的重叠
  # spaa::niche.overlap(datasample[,1:3], method = "levins")
  # #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247486549&idx=1&sn=05c283e4cf0afc4d338f1ca1c11dc5fc&chksm=ea71f29edd067b8829272ee565242ad9b08d18c4813d9e025144b41b9506c4cba40b754e2b6e&scene=21#wechat_redirect
  # lib_ps("MicroNiche",library = FALSE)
  # ?MicroNiche::levins.Bn()
  # ?MicroNiche::hurlberts.Bn()
  # MicroNiche::hurlberts.Bn(otutab)
  # ?MicroNiche::feinsingers.PS()
  # #https://mp.weixin.qq.com/s?__biz=MzI1OTk3NjEwMA==&mid=2247485087&idx=1&sn=95d121fc4080aa9af2a61de22694093d&chksm=ea71f854dd067142c0b201657c57aeb6102cef5f738fc89380783efa8c413580515b8536e5a6&scene=21#wechat_redirect
  # lib_ps("indicspecies",library = FALSE)
  #
  # #labdsv包的indval函数可计算群落中的指示物种
  # lib_ps("labdsv",library = FALSE)
  # labdsv::indval()
  # #分析微生物群落和环境因子相关性的工具vegan::bioenv
  # vegan::bioenv()
  # #表征微生物抵抗性的一个指标(要有实验-对照的关系才能算)
  # #Rs=1-(2*|D0|/(C0+|D0|))
  # #其中C0为对照组样本的α多样性指数，D0为对照组样本多样性和实验组样本α多样性之间的差值。
  # #Rs取值为-1到1。越接近1表明处理效应非常小，即微生物抵抗性很强；Rs越小表明处理效应越强，即微生物抵抗力越弱。
  #
  # #DNCI群落构建
  # lib_ps("DNCImper",library = FALSE)
  # ?PerSIMPER
  # #用法
  # DNCImper:::PerSIMPER(
  #   matrixSIMP,  #行为样本，列为物种
  #   Groups,      #只允许两个组
  #   leg = FALSE, #是否添加图例
  #   count = TRUE,#显示完成置换的次数
  #   dataTYPE = "prab", #prab为0-1数据；count为丰度数据
  #   Nperm = 1000,#置换次数
  #   plotSIMPER = TRUE
  # )
  # A <- DNCImper:::PerSIMPER(Matrix, Group,Nperm = 100,count = FALSE)
  # #上述只能针对两个组，若有更多的组，可用PerSIMPER_overall，计算整体的PerSIMPER
  # #DNCI.ese  计算DNCI效应量,只允许两组。
  # #三组或以上用DNCI.ses_overall，计算整体的DNCI。
  # #三组或以上若用DNCI_multigroup，计算两两成对DNCI。
}


#======Life Game======

#' Life Game Simulation
#'
#' @param file gif filename
#' @param time how many times the life game continue.
#'
#' @return a gif file
#' @export
life_game=function(file="LifeGame",time = 100){
  # Game of Life
  # Refer to: https://zhuanlan.zhihu.com/p/136727731
  lib_ps("animation",library = FALSE)
  lib_ps("pheatmap",library = FALSE)   #加载pheatmap

  ### 构造初始状态：
  set.seed(2022-2-21)
  size = 20  # 矩阵的行和列数
  d = round(runif(size*size,0,0.6)) # 最大值低一些，保证初始有值的少一些。
  start = matrix(data=d,ncol=size,nrow=size)

  ### 求下一个状态时格子周围值的和：
  fun = function(m,x,y,size){ # m为当前状态的矩阵；x和y为坐标；size为矩阵大小
    fun.sum = 0
    for (i in c(x-1,x,x+1)){ #依次遍历一个格子周围3x3的邻居格子
      for (j in c(y-1,y,y+1)){
        #如果格子在角落或者边，则邻居的值直接为0
        if (i  > 0 & i <= size & j > 0 & j <= size)  fun.sum = fun.sum + m[i,j]  # 把9个格子先求和
      }
    }
    fun.sum = fun.sum - m[x,y]  # 减去中间格子的值，即为周围8个值的和
  }
  ### 设置运行次数
  time = time
  life = list()
  life[[1]] = start
  for (k in 2:time){ # k = 3
    life.next = matrix(data = 0,ncol=size,nrow = size)
    for (i in 1:size){
      for (j in 1:size){
        fun.sum = fun(life[[k-1]],i,j,size) # 运行上述函数

        # 判断下个状态时当前位置是否有值存在。判断依据来自 http://nonoas.gitee.io/webproj/LifeGame/

        #孤单死亡：如果细胞的邻居小于等于1个，则该细胞在下一次状态将死亡；

        #拥挤死亡：如果细胞的邻居在4个及以上，则该细胞在下一次状态将死亡；

        #稳定：如果细胞的邻居为2个或3个，则下一次状态为稳定存活；

        #复活：如果某位置原无细胞存活，而该位置的邻居为2个或3个，则该位置将复活一个细胞。

        life.next[i,j] =  ifelse( fun.sum == 2|fun.sum == 3, 1, 0)
      }

    }
    life[[k]] = life.next
  }
  #输出gif：
  animation::saveGIF(
    (for (k in 1:time){
      pheatmap::pheatmap(life[[k]],cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=FALSE,
               legend = FALSE,cellwidth = 20, cellheight = 20,
               color = c("white","black"),main = "Life Game")
    }), movie.name = paste0(file,".gif")
  )
}


