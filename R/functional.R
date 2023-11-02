#=====functional gene analysis====
#
#' Gene symbolid transfer to entrezIDs (human gene)
#'
#' @param genes gene symbols e.g:ASGR2
#'
#' @return gene entrezIDs dataframe
#' @export
#'
#' @examples
#' genes=c("ASGR2","BEST1","SIGLEC16","ECRP","C1QC","TCN2","RNASE2",
#'     "DYSF","C1QB","FAM20A","FCGR1A","CR1","HP","VSIG4","EGR1")
#' gene2id(genes)->geneid
gene2id<-function(genes){
  lib_ps("clusterProfiler","org.Hs.eg.db",library = F)
  #entrezIDs=AnnotationDbi::mget(genes,org.Hs.egSYMBOL2EG, ifnotfound=NA)  # 找出基因对应的ID
  suppressMessages({entrezIDs=AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENTREZID")})
  entrezIDs=as.character(entrezIDs) # 获取数据
  rt=data.frame(genes,entrezID=entrezIDs)       # 添加一列entrezID
  #rt=rt[(rt[,"entrezID"])!="NA",]   # 删除没有基因的ID
  rt=rt[!is.na(rt[,"entrezID"]),]   # 删除没有基因的ID
  return(rt)
}

#' prepare the GEO data
#'
#' @param my_id GEOid
#' @param GEO_dir GEO download dir
#' @param file the downloaded file
#'
#' @export
#'
pre_GEO=function(my_id,GEO_dir="GEO_data",file=NULL){
  #1.Importing the data
  ## change my_id to be the dataset that you want.
  lib_ps("GEOquery","Biobase",library = F)

  if(is.null(file)){
    if (file.exists(paste0(GEO_dir,"/",my_id,"_series_matrix.txt.gz"))) {
      file=paste0(GEO_dir,"/",my_id,"_series_matrix.txt.gz")
    }
  }

  if(is.null(file)){
    gse <- GEOquery::getGEO(my_id[1],destdir = GEO_dir)
    gse<-gse[[1]]
  }
  else gse <- GEOquery::getGEO(filename = file, getGPL = F)

  GPL_version= gse@annotation

  meta=Biobase::pData(gse) ## print the sample information
  GPL_data_11=Biobase::fData(gse) ## print the gene annotation
  GSE_data_expr=Biobase::exprs(gse) ## print the expression data

  # 进行注释------
  GSE_data_expr <- GSE_data_expr %>% as.data.frame()%>%
    tibble::rownames_to_column() %>%
    reshape::rename(c(rowname = "prode_id"))%>%mutate(prode_id=as.character(prode_id))

  if(GPL_version%in%c("GPL570","GPL571","GPL96","GSE22873","GPL1261","GPL81")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,`Gene Symbol`) %>%
      reshape::rename(c(ID = "prode_id",`Gene Symbol` = "GeneSymbol"))
    GPL$GeneSymbol=gsub(" ///.*","",GPL$GeneSymbol) #一对多取第一个？
  }
  else if(GPL_version%in%c("GPL6254")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,`GENE_SYMBOL`) %>%
      reshape::rename(c(ID = "prode_id",`GENE_SYMBOL` = "GeneSymbol"))
    GPL$GeneSymbol=gsub(" ///.*","",GPL$GeneSymbol) #一对多取第一个？
  }
  else if(GPL_version%in%c("GPL10787")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,`GENE_SYMBOL`) %>%
      reshape::rename(c(ID = "prode_id",`GENE_SYMBOL` = "GeneSymbol"))
  }
  else if (GPL_version%in%c("GPL6947","GPL6883","GPL4866","GPL4006")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,`Symbol`) %>%
      reshape::rename(c(ID = "prode_id",`Symbol` = "GeneSymbol"))
    GPL$GeneSymbol=gsub("/.*","",GPL$GeneSymbol) #一对多取第一个？
  }
  else if (GPL_version%in%c("GPL6244","GPL17586","GPL6246")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,gene_assignment) %>%
      reshape::rename(c(ID = "prode_id",gene_assignment = "GeneSymbol"))
    tmp<-strsplit(GPL$GeneSymbol,split = " // ")
    lapply(tmp,\(x)x[2])%>%do.call(c,.)->GPL$GeneSymbol
    na.omit(GPL)->GPL
  }
  else if (GPL_version%in%c(GPL1536)){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,GENE) %>%
      reshape::rename(c(ID = "prode_id",GENE = "GeneSymbol"))
    GPL$GeneSymbol=gsub("/.*","",GPL$GeneSymbol)
  }
  else if (GPL_version%in%c("GPL10739")){
    GPL <- GPL_data_11 %>%
      dplyr::select(ID,gene_assignment) %>%
      reshape::rename(c(ID = "prode_id",gene_assignment = "GeneSymbol"))
    tmp<-strsplit(GPL$GeneSymbol,split = " // ")
    lapply(tmp,\(x)x[2])%>%do.call(c,.)->GPL$GeneSymbol
    na.omit(GPL)->GPL
    tmp<-strsplit(GPL$GeneSymbol,split = " /// ")
    lapply(tmp,\(x)x[1])%>%do.call(c,.)->GPL$GeneSymbol
    na.omit(GPL)->GPL
  }
  else stop("unknown GPL version, please check and modify the code.")

  GPL <-GPL%>%
    add_count(GeneSymbol,name = "n_gene") %>%
    add_count(prode_id,name = "n_prode_id") %>%
    dplyr::filter(n_gene < 100) %>% # 选择能有对应基因的探针,太多了说明非特异
    dplyr::select(c(-3,-4))%>%mutate_all(as.character)

  GSE_data_expr[is.na(GSE_data_expr)]=0
  GSE_expr <- GSE_data_expr %>%
    inner_join(GPL, ., by = "prode_id") %>% #p2s和上一步得到的结果再取交集，p2s放在右边是以它为准
    dplyr::select(-1) %>%  #去除第一列probe_id
    data.frame() %>% #因为aggregate需要数据框格式
    aggregate(.~ GeneSymbol, data = ., mean) %>%  #以symbol作为因子水平，把相似的数据放在一起取均值，最大值max，中位值median
    tibble::column_to_rownames(var = "GeneSymbol")

  rownames(GSE_expr)%>%gsub(" ","",.)->rownames(GSE_expr)

  saveRDS(list(meta=meta,GSE_expr=GSE_expr),file = paste0(GEO_dir,"/",my_id,".RDS"))
  return(list(meta=meta,GSE_expr=GSE_expr))
}

