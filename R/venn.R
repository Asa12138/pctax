#venn========
#' plot a general venn (upset, flower)
#'
#' @param object list/data.frame/pc_otu
#' @param ... additional
#'
#' @return a plot
#' @export
#'
#' @examples
#'aa=list(a=1:3,b=3:7,c=2:4)
#'venn(aa)
venn <- function(object,...){
  UseMethod("venn", object)
}

#' Calculate dataframe rownames(>0) in each columns
#'
#' @param otu_time a dataframe
#'
#' @return a list
#'
#' @examples
#'data.frame(row.names = letters[1:10],col1=sample(0:1,10,T),col2=sample(0:1,10,T))->aa
#'venn_cal(aa)
#' @export
venn_cal<-function(otu_time){
  aa=list()
  for (i in 1:ncol(otu_time)){
    name=colnames(otu_time)[i]
    aa[[name]]=rownames(otu_time[otu_time[,i]>0,])
  }
  return(aa)
}

#' @method venn list
#' @rdname venn
#' @param mode "venn","venn2","upset","flower"
#' @export
venn.list<-function(aa,mode="venn",...){
  if(length(aa)>4&&mode=="venn")print("venn < 4, recommend upset or flower")
  if(mode=="venn")lib_ps("ggvenn");ggvenn::ggvenn(aa)->p
  if(mode=="venn2"){
    lib_ps("Vennerable")
    Vennerable::Venn(aa)->aap
    plot(aap,...)
    # plot(aap,type="triangles")
    # plot(aap, doWeights = FALSE)
    # plot(aap, doWeights = FALSE,type="ellipses")
    # plot(aap, doWeights = FALSE,type="ChowRuskey")
  }
  if(mode=="upset"){
    lib_ps("UpSetR")
    UpSetR::upset(UpSetR::fromList(aa), order.by = "freq",nsets = length(aa),nintersects = 30)->p
    }
  if(mode=="flower"){
    lib_ps("RColorBrewer","plotrix")
    otu_num=length(aa[[1]])
    core_otu_id=aa[[1]]
    for (i in 2:length(aa)){
      core_otu_id <- intersect(core_otu_id, aa[[i]])
      otu_num <- c(otu_num, length(aa[[i]]))
    }
    core_num <- length(core_otu_id)
    otu_num<-otu_num-core_num
    sample_id<-names(aa)
    n <- length(sample_id)

    ellipse_col <- colorRampPalette(brewer.pal(10,"Set3"))(n)#椭圆的颜色设置
    start = 90; a = 0.5; b = 2.2; r = 0.5; ellipse_col = ellipse_col; circle_col = 'white'

    par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
    #start为椭圆的起始位置, a为椭圆的短轴, b为椭圆的长轴, r为中心圆的半径
    plot(c(0,10),c(0,10),type='n')
    deg <- 360 / n
    res <- lapply(1:n, function(t){
      draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), #绘制椭圆
                   y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                   col = ellipse_col[t],
                   border = ellipse_col[t],
                   a = 0.6, b = 2.2, angle = deg * (t - 1))

      text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),#在图上设置每个样品的otu_num数目
           y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
           otu_num[t])

      if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),#设置每个样品名
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample_id[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = 1
        )
      } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample_id[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = 1
        )
      }
    })
    draw.circle(x = 5, y = 5, r = 1.3, col = circle_col, border = NA)#绘制中心圆
    text(x = 5, y = 5, paste('Core:', core_num))#设置中心圆名字
    #sample_id = sample_id, otu_num = otu_num, core_num = core_num,
    #start = 90, a = 0.5, b = 2.2, r = 0.5, ellipse_col = ellipse_col, circle_col = 'white'
    p=NULL
  }
  return(p)
}

#' @method venn data.frame
#' @rdname venn
#' @export
venn.data.frame<-function(otutab,mode="venn"){
  venn_cal(otutab)->aa
  venn.list(aa,mode = mode)
}

#' @method venn pc_otu
#' @rdname venn
#' @export
venn.pc_otu<-function(pc,mode="venn",tbl="otutab"){
  pc_valid(pc)
  otutab<-pc$tbls[[tbl]]
  venn.data.frame(otutab,mode = mode)
}
