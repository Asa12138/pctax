#data trans========
#数据转化、筛选方法

#' Transfer your data
#' @param df dataframe
#' @param method "cpm","minmax","asinh","total", "max", "frequency", "normalize", "range", "rank", "rrank",
#' "standardize", "pa", "chi.square", "hellinger", "log", "clr", "rclr", "alr"
#' @param margin 1 for row and 2 for column(default: 2)
#' @param ... additional
#' @examples
#'data(otutab)
#'trans(otutab,method="cpm")
#' @seealso \code{\link[vegan]{decostand}}
trans<-function(df,method = "normalize",margin=2,...){
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

#' Filter your data
#'
#' @param tab dataframe
#' @param sum the rowsum should bigger than sum(default:10)
#' @param exist the exist number bigger than exist(default:1)
#'
#' @return input object
#' @export
#'
#' @examples
#'data(otutab)
#'guolv(otutab)
guolv<-function(tab,sum=10,exist=1){
  tab[rowSums(tab)>sum,]->tab
  tab[rowSums(tab>0)>exist,]->tab
  return(tab)
}

#' Group your data
#'
#' @param otutab dataframe
#' @param group group vector
#' @param margin 1 for row and 2 for column(default: 2)
#' @param act do (default: mean)
#' @rdname hebing
#' @return input object
#' @export
#'
#' @examples
#' data(otutab)
#' hebing(otutab,metadata$Group)
hebing<-function(otutab,group,margin=2,act='mean'){
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

