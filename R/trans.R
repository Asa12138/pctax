#data trans========
#数据转化、筛选方法

#' Transfer your data
#'
#' @param object object
#' @param ... additional
#'
#' @return input object
#' @export
#'
#' @examples
#'data(otutab)
#'trans(otutab,method="cpm")
trans <- function(object,...){
  UseMethod("trans", object)
}


#' @param df dataframe
#' @param method "cpm","minmax","asinh","total", "max", "frequency", "normalize", "range", "rank", "rrank",
#' "standardize", "pa", "chi.square", "hellinger", "log", "clr", "rclr", "alr"
#' @param margin 1 for row and 2 for column(default: 2)
#' @param ... additional
#' @rdname trans
#' @exportS3Method
#' @seealso \code{\link[vegan]{decostand}}

trans.data.frame<-function(df,method = "normalize",margin=2,...){
  all=c("cpm","minmax","acpm","total", "max", "frequency", "normalize", "range", "rank", "rrank",
        "standardize", "pa", "chi.square", "hellinger", "log", "clr", "rclr", "alr")
  if (!method%in%all)stop("methods should be one of ",all)
  if(method=="cpm"){
    df=apply(df, margin, \(x){x*10**6/sum(x)})
  }
  else if(method=="minmax"){
    df=apply(df,margin,mmscale,...)
  }
  else if(method=="acpm"){
    df=asinh(apply(df, margin, \(x){x*10**6/sum(x)}))
  }
  else df=vegan::decostand(df,method = method,margin,...)
  return(data.frame(df,check.names = F))
}

#' @rdname trans
#' @exportS3Method
trans.pc_otu<-function(pc,method = "normalize",tbl="otutab",margin=2,...){
  pc_valid(pc)
  df<-pc$tbls[[tbl]]
  pc$tbls[[method]]=trans.data.frame(df,method = method,margin=margin,...)
  return(pc)
}

#' Filter your data
#'
#' @param object object
#' @param ... additional
#'
#' @return input object
#' @export
#'
#' @examples
#'data(otutab)
#'guolv(otutab)
guolv <- function(object,...){
  UseMethod("guolv", object)
}

#' @param tab dataframe
#' @param sum the rowsum should bigger than sum(default:10)
#' @param exist the exist number bigger than exist(default:1)
#'
#' @rdname guolv
#' @exportS3Method
guolv.data.frame<-function(tab,sum=10,exist=1){
  tab[rowSums(tab)>sum,]->tab
  tab[rowSums(tab>0)>exist,]->tab
  return(tab)
}

#' @rdname guolv
#' @exportS3Method
guolv.pc_otu<-function(pc,sum=10,exist=1){
  pc_valid(pc)
  df<-pc$tbls$otutab
  pc$tbls[["filter"]]=guolv.data.frame(df,sum=sum,exist=exist)
  return(pc)
}

#' Group your data
#'
#' @param object object
#' @param ... additional
#'
#' @return input object
#' @export
#'
#' @examples
#' data(otutab)
#' hebing(otutab,metadata$Group)
hebing <- function(object,group,...){
  UseMethod("hebing", object)
}

#' @param otutab dataframe
#' @param group group vector
#' @param margin 1 for row and 2 for column(default: 2)
#' @param act do (default: mean)
#' @rdname hebing
#' @exportS3Method
hebing.data.frame<-function(otutab,group,margin=2,act='mean'){
  if (margin==2) {
    aggregate(t(otutab),FUN=act,by=list(factor(group)))->a
    a[,-1]->a
    data.frame(t(a))->a
    levels(factor(group))->colnames(a)
  }
  else{
    aggregate(otutab,FUN=act,by=list(factor(group)))->a
    a[,-1]->a
    levels(factor(group))->rownames(a)
  }
  return(a)
}

#' @rdname hebing
#' @exportS3Method
hebing.pc_otu<-function(pc,group,margin=2,act='mean',tbl="otutab"){
  pc_valid(pc)
  otutab<-pc$tbls[[tbl]]
  if(!group%in%colnames(pc$metas$metadata))stop(group,"not found in metadata")
  group1=pc$metas$metadata[[group]]
  pc$tbls[[group]]=hebing.data.frame(otutab,group=group1,margin=margin,act=act)
  return(pc)
}

