#rare==========

#' Rarefy a otutab
#'
#' @param otutab otutab
#' @param sample number
#'
#' @return a rarefied otutab
#' @export
#'
#' @examples
#' data(otutab)
#' rarefaction(otutab)
rarefaction<-function(otutab,sample=NULL){
  lib_ps("vegan")
  if(is.null(sample))sample=min(colSums(otutab))
  otutab1=rrarefy(t(otutab),sample = sample)
  return(as.data.frame(t(otutab1)))
}

#' Rare the sample and plot
#'
#' @param otutab otutab
#' @param rep repeats number
#' @param count_cutoff cutoff to be 0
#'
#' @export
#'
#' @examples
#' rare_curve_sam(otutab)
rare_curve_sam<-function(otutab, rep=30, count_cutoff=1){
  # 默认值绘制样本箱线图稀释曲线
  count = otutab
  count[count < count_cutoff] = 0
  count[count >= count_cutoff] = 1
  n = dim(count)[2]
  x = unique(as.integer(seq(1, n, length.out = n)))
  result = data.frame(sample = c(0), richness = c(0))
  for (i in x) {
    m = choose(n, i)
    if (m > rep) {
      m = rep
    }
    for (j in 1:m) {
      idx = sample(n, i)
      temp = count[, idx, drop = F]
      temp1 = temp[rowSums(temp) > 0, , drop = F]
      result = rbind(result, c(i, dim(temp1)[1]))
    }
  }
  result = result[-1, ]
  result$sample = factor(result$sample, levels = unique(result$sample))
  p = ggplot(result, aes(x = sample, y = richness, fill = sample)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.5,size = 0.1) + xlab("Sample size") +
    ylab("Feature richness") + theme(legend.position = "none",
                                     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  return(list(result,p))
}

#' Plot a rare curve
#'
#' @param otutab otutab
#' @param step default 20000
#' @param method one of "richness","chao1","ace","gc","shannon","simpson","pd","pielou"
#' @param mode 1 for little table, 2 for big
#' @param threads use how many threads to calculate (default:4)
#' @export
#'
#' @examples
#' data(otutab)
#' rare_curves(otutab,mode=1)
rare_curves<-function(otutab,step = 40000,method = 'richness',mode=2,threads=4){
  #根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
  alpha_curves <- \(x, step, method = 'richness', tree = NULL) {
    library(pcutils)
    x_nrow <- nrow(x)
    rare <- rowSums(x)

    alpha_rare <- data.frame()
    for (i in 1:x_nrow) {
      step_num <- seq(0, rare[i], step)
      if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
      alpha_rare_i <- NULL
      for (step_num_n in step_num){
        pctax::a_diversity(data.frame(t(rrarefy(x[i,],step_num_n))), method = method, tree = tree)->tmp_i
        alpha_rare_i <- c(alpha_rare_i,unlist(tmp_i))
      }
      alpha_rare <- rbind(alpha_rare,
                          data.frame(rare=step_num,alpha=alpha_rare_i,sample=rownames(x)[i]))
    }
    alpha_rare
  }

  if (mode==1){
    t(otutab)->otu
    step = 2000
  }else t(round(trans(otutab,"cpm"),0))->otu

  lib_ps("parallel","doSNOW")
  pb <- txtProgressBar(max =5, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  plot_richness1 <- foreach(randomize = 1:5,.options.snow = opts,.packages = c("vegan")) %dopar% {
    alpha_curves(otu, step = step, method = 'richness')
  }
  stopCluster(cl)
  plot_richness<-do.call(rbind,plot_richness1)

  df<-plot_richness%>%group_by(.,rare,sample)%>%summarise(.,mean(alpha),sd(alpha))%>%dplyr::select(everything(),alpha=3,sd=4)
  df$sample<-df$sample%>%factor(.,levels=colnames(otutab))
  nn=nlevels(df$sample)

  p=ggplot(df, aes(x=as.numeric(rare), y=alpha, color = sample)) +geom_line()+
    geom_errorbar(aes(ymin=alpha-sd,ymax=alpha+sd))+
    geom_point(size=1)+labs(x=ifelse(mode==1,"reads(count)",'reads(cpm)'),y=method)+theme_bw()

  return(list(df=df,p=p))
}

#'Calculate Abundance-occupancy_relationship
#'
#'@param object object, such as data.frame or pc_otu
#'@export
#'@return AOR
#'@examples
#'data("otutab")
#'aor(otutab)->AOR
#'plot(AOR)

aor <- function(object,...){
  UseMethod("aor")
}

#' @param otutab otutab
#'
#' @param top_r percentage of top relative abundance
#' @param ocup_n percentage of top occupied
#'
#' @method aor data.frame
#' @import dplyr
#' @rdname aor
#' @exportS3Method
aor.data.frame<-function(otutab,top_r=0.7,ocup_n=ceiling(0.8*ncol(otutab))){
  lib_ps("dplyr")
  otutab%>%transmute(abundance=rowSums(.),relative_abund=abundance/sum(abundance),
                     occupancy=rowSums(.>0))%>%arrange(-abundance)->AOR
  AOR%>%mutate(log_abund=log10(relative_abund))->AOR

  y=AOR$occupancy;x=AOR$log_abund
  tmp<-try({
    fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y))
    #summary(fit)
    AOR$fit_y=predict(fit, newdata = data.frame(x = x))
  })
  if("try-error" %in% class(tmp)){
    AOR$fit_y=0
  }
  #find core otus
  AOR%>%mutate(cum=cumsum(relative_abund),core=ifelse((cum<top_r)&(occupancy>ocup_n),"core","others"))->AOR
  class(AOR)=c("AOR","data.frame")
  return(AOR)
}

#' @method aor pc_otu
#' @rdname aor
#' @exportS3Method
aor.pc_otu<-function(pc,top_r=0.7,ocup_n=0.8,tbl="otutab"){
  pc_valid(pc)
  otutab<-pc$tbls[[tbl]]
  pc$otus$AOR<-aor.data.frame(otutab,top_r=top_r,ocup_n=ceiling(ocup_n*ncol(otutab)))
  return(pc)
}


#' Plot a AOR
#' @param AOR a AOR object
#' @import scales
#' @rdname aor
#' @exportS3Method
plot.AOR<-function(AOR){
  lib_ps("scales")
  p1<-ggplot()+
    geom_point(data = AOR,aes(x=10^log_abund,y=occupancy),col="#79A5C9",alpha=0.6,size=2)+
    geom_line(data = AOR,aes(x=10^log_abund,y=fit_y),col="red",size=1.5)+
    scale_x_log10(labels =scales::label_log())+
    theme_bw()+
    labs(x="mean relative abundnace (%)",y="occupancy")
  #scale_x_continuous()

  p2<-ggplot()+geom_point(data = AOR,aes(x=10^log_abund,y=occupancy,col=core),alpha=0.6,size=2)+
    scale_x_log10(labels =scales::label_log())+
    scale_color_manual(values = c("others"="grey","core"="#79A5C9"))+
    theme_bw()+
    labs(x="mean relative abundnace (%)",y="occupancy")
  return(list(p1,p2))
}
