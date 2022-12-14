#a_diversity==========

#' Calculate a_diversity of otutab
#' @export
#'
#' @examples
#' data(otutab)
#' a_diversity(otutab)->a_res
#' plot(a_res,metadata,"Group")
#' plot(a_res,metadata,"env1")
a_diversity <- function(object,...){
  UseMethod("a_diversity", object)
}

#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param method one of "all","richness","chao1","ace","gc","shannon","simpson","pd","pielou"
#' @param tree a iphylo object match the rownames of otutab
#' @param digits maintance how many digits
#' @rdname a_diversity
#' @return a a_res object
#' @exportS3Method
a_diversity.data.frame<-function(otutab,method=c("simpson","shannon"),tree=NULL,digits=4){
  lib_ps("vegan")
  all=c("all","richness","chao1","ace","gc","shannon","simpson","pd","pielou")
  if(!all(method%in%all))stop(paste0("methods should be some of ",paste0(all,collapse = ",")))
  if("all"%in%method)method=all[-1]
  x=t(otutab)
  a_res=data.frame(row.names = colnames(otutab))
  if("richness"%in%method){Richness <-rowSums(x>0);a_res=cbind(a_res,Richness)}
  if("chao1"%in%method){Chao1 <- vegan::estimateR(x)[2, ];a_res=cbind(a_res,Chao1)}
  if("ace"%in%method){ACE <- vegan::estimateR(x)[4, ];a_res=cbind(a_res,ACE)}
  if("gc"%in%method){Goods_Coverage <- 1 - rowSums(x <= 1) / rowSums(x);a_res=cbind(a_res,Goods_Coverage)}
  if("shannon"%in%method){Shannon <- vegan::diversity(x, index = 'shannon');a_res=cbind(a_res,Shannon)}
  #注意，这里是Gini-Simpson 指数
  if("simpson"%in%method){Simpson <- vegan::diversity(x, index = 'simpson');a_res=cbind(a_res,Simpson)}
  if("pielou"%in%method){
    Pielou_evenness<-vegan::diversity(x, index = 'shannon')/log(rowSums(x>0))
    a_res=cbind(a_res,Pielou_evenness)
  }
  if("pd"%in%method){
    if(is.null(tree))warning("pd need tree!")
    else{
      lib_ps("picante")
      match.phylo.comm(tree,x)->match_p
      pds <- picante::pd(match_p$comm, match_p$phy, include.root = FALSE)
      PD <- pds[ ,1]
      a_res=cbind(a_res,PD)
      #净相关指数
      # NRI=-ses.mpd(x,cophenetic(spe_nwk),null.model="taxa.labels")[6]
      # names(NRI) <- 'NRI'
      #最近邻体指数
      # NTI=-ses.mntd(x,cophenetic(spe_nwk),null.model="taxa.labels")[6]
      # names(NTI) <- 'NTI'
      # result <- cbind(result, PD_whole_tree,NRI,NTI)
    }
  }
  a_res<-round(a_res,digits)
  class(a_res)<-c("a_res","data.frame")
  return(a_res)
}

#' @exportS3Method
#'
#' @rdname a_diversity
a_diversity.pc_otu<-function(pc,method="all",tbl="otutab"){
  pc_valid(pc)
  otutab<-pc$tbls[[tbl]]
  pc$metas$a_res<-a_diversity.data.frame(otutab,method = method)
  return(pc)
}

#' @exportS3Method
#'
#' @rdname a_diversity
a_diversity.numeric<-function(x,method=c("simpson","shannon"),tree=NULL){
  return(a_diversity(data.frame(Sample=x)))
}

#' Plot a_res object
#'
#' @param a_res a a_res object
#' @param metadata metadata
#' @param group one of colname of metadata
#' @param ... addditional parameters for \code{\link{group_box}} or \code{\link{my_lm}}
#'
#' @return patchwork object,you can change theme with &
#' @export
#'
#' @seealso \code{\link{a_diversity}}
#'
plot.a_res<-function(a_res,metadata,group,...){
  a_res<-a_res[rownames(metadata),,drop=F]
  group1=metadata[,group]
  if(is.numeric(group1)&!is.factor(group1)){
    p=my_lm(a_res,group,metadata,...)
  }
  else {
    p=group_box(a_res,group,metadata,...)
    }
  return(p)
}


#test phylogenetic diversity
if(F){
  lib_ps("picante")
  data("phylocom")
  View(phylocom$sample)
  ggtree(phylocom$phylo)+geom_tiplab()+theme_tree2()
  samp=data.frame(a=1:3,b=2:4,c=c(1,2,0),d=c(0,0,3),e=0,row.names = paste0("plot",1:3))
  read.tree(text = "(c:2,(a:1,b:1):2,d:1,f:1):1;",)->test
  ggtree(test)+geom_tiplab()+theme_tree2()
  #prune.sample(samp,test)->test_prune
  match.phylo.comm(test,samp)->match_p
  match_p$phy->test;match_p$comm->samp

  test_prune%>%ggtree(.)+geom_tiplab()+theme_tree2()
  pd(samp,test,include.root = F)
  pd(samp,test_prune,include.root = F)

  #picante::ses.pd(samp,test)
  #tree noedes distance
  cophenetic(test)
  #每个样方有的物种对的平均谱系距离mpd
  mpd(samp,cophenetic(test))
  #随机化mpd
  mpd(samp,taxaShuffle(cophenetic(test)))
  #ses.mpd直接计算了
  ses.mpd(samp,cophenetic(test))->mpd_ses
  #净相关指数nri,>0聚集，<0发散
  mpd_ses%>%mutate(nri=-1*(mpd.obs-mpd.rand.mean)/mpd.rand.sd)%>%select(nri)
  #nti类似nri
  #mnpd最近谱系距离均值
  mntd(samp,cophenetic(test))
  ses.mntd(samp,cophenetic(test))->mnpd_ses
  mnpd_ses%>%mutate(nti=-1*(mntd.obs-mntd.rand.mean)/mntd.rand.sd)%>%select(nti)

  #beta-mpd
  comdist(samp,cophenetic(test))
  #beta-mntd
  comdistnt(samp,cophenetic(test))
  #pcd
  picante::pcd(samp,test)
  #phylosor
  picante::phylosor(samp,test)
  picante::psd(samp,test)
  picante::raoD(samp,test)
  picante::unifrac(samp,test)

}
