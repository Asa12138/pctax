#'Create a pc_otu class object
#'
#'@param otutab an otutab data.frame, samples are columns, taxs are rows.
#'@param metadata a metadata data.frame, samples are rows
#'@param taxonomy a taxomomy data.frame, look out the rowname of taxonomy and otutab should matched!
#'@export
#'@return pc_otu
#'@examples
#'data("otutab")
#'otu1<-pc_otu(otutab,metadata)
#'otu1
pc_otu<-function(otutab=data.frame(),metadata=data.frame(),taxonomy=NULL,...){
  if(!is.data.frame(otutab))stop("otutab must be a df!")
  if(!is.data.frame(metadata))stop("metadata must be a df!")
  if(!is.null(taxonomy)){
    if(!is.data.frame(taxonomy))stop("taxonomy must be a df!")
    if(length(taxonomy)!=7)stop("taxonomy should be a df with 'k,p,c,o,f,g,s'!")
    fillNAtax(taxonomy)->taxonomy
    taxonomy[match(rownames(otutab),rownames(taxonomy)),]%>%na.omit()->taxonomy
  }
  ids=intersect(rownames(metadata),colnames(otutab))
  if(length(ids)<1)stop("no samples left, check the rownames of metadata or colnames of otutable")
  metadata=metadata[ids,,drop=F]
  otutab=otutab[,ids,drop=F]
  otutab=otutab[rowSums(otutab)>0,,drop=F]
  pc<-list(
    tbls=list(otutab=otutab),
    metas=list(metadata=metadata),
    otus=list(taxonomy=taxonomy),...
  )
  class(pc)<-append("pc_otu",class(pc))
  return(pc)
}

#' Judge pc_otu is valid or not
#' @param pc a pc_otu object
#' @export
pc_valid<-function(pc){
  for (i in pc$tbls)if(!is.data.frame(i) & !is.null(i))stop("tbls must be df!")
  for (i in pc$metas)if(!is.data.frame(i) & !is.null(i))stop("metas must be df!")
  for (i in pc$otus)if(!is.data.frame(i) & !is.null(i))stop("otus must be df!")
  if(!all(rownames(pc$metas$metadata)==colnames(pc$tbls$otutab)))stop("Wrong samples! check the rownames of metadata or colnames of otutable")
  # print("it's OK!")
  return(T)
}

#'@method print pc_otu
#'@export
print.pc_otu<-function(pc){
  sprintf("There are %d otus and %d samples!",nrow(pc$tbls$otutab),ncol(pc$tbls$otutab))%>%print()
  for (i in names(pc)){
    dabiao(i)
    if(i%in%c("tbls","metas","otus")){
      for (j in names(pc[[i]])){
        dabiao(j,40)
        if("data.frame"%in%class(pc[[i]][[j]]))print(head(pc[[i]][[j]]))
        else print(pc[[i]][[j]])
      }
    }
    else {
      if("data.frame"%in%class(pc[[i]]))print(head(pc[[i]]))
      else print(pc[[i]])}
  }
}

#' @method summary pc_otu
#'@export
summary.pc_otu<-function(pc){
  pc_valid(pc)
  sprintf("There are %d otus and %d samples!",nrow(pc$tbls$otutab),ncol(pc$tbls$otutab))%>%print()
  dabiao("tables, include some data filter and tranformat")
  print(names(pc$tbls))
  dabiao("metadatas, include some statistics")
  print(names(pc$metas))
  dabiao("otus annotation, include some statistics")
  print(names(pc$otus))
  dabiao("some other indexs:")
  print(names(pc)[-1:-3])
}

#'@title Add taxonomy for a pc_otu object
#'@export
#'@param pc a pc_otu object
#'@param taxonomy a taxomomy data.frame, look out the rownames of taxonomy and otutab should matched!
#'@examples
#'data("otutab")
#'otu1<-pc_otu(otutab,metadata)
#'otu1<-add_tax(otu1,taonomy)
add_tax<-function(pc,taxonomy){
  if(!"pc_otu"%in%class(pc))stop(pc,"should be a pc_otu")
  if(!is.data.frame(taxonomy))stop("taxonomy must be a df!")
  if(length(taxonomy)!=7)stop("taxonomy should be a df with 'k,p,c,o,f,g,s'!")

  fillNAtax(taxonomy)->taxonomy
  taxonomy[match(rownames(pc$tbls$otutab),rownames(taxonomy)),]%>%na.omit()->taxonomy
  pc$otus$taxonomy=taxonomy
  return(pc)
}

