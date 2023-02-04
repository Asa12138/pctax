#=====functional gene analysis
#

#' Gene symbol transfer to entrezIDs
#'
#' @param genes gene symbols e.g:ASGR2
#'
#' @return gene entrezIDs dataframe
#' @export
#'
#' @examples
#' genes=c("ASGR2","BEST1","SIGLEC16","ECRP","C1QC","TCN2","RNASE2","DYSF","C1QB","FAM20A","FCGR1A","CR1","HP","VSIG4","EGR1")
#' gene2id(genes)->geneid
gene2id<-function(genes){
  lib_ps("clusterProfiler","org.Hs.eg.db")
  #entrezIDs=AnnotationDbi::mget(genes,org.Hs.egSYMBOL2EG, ifnotfound=NA)  # 找出基因对应的ID
  entrezIDs=mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENTREZID")
  entrezIDs=as.character(entrezIDs) # 获取数据
  rt=data.frame(genes,entrezID=entrezIDs)       # 添加一列entrezID
  #rt=rt[(rt[,"entrezID"])!="NA",]   # 删除没有基因的ID
  rt=rt[!is.na(rt[,"entrezID"]),]   # 删除没有基因的ID
  return(rt)
}

#' GO enrich or KEGG enrich
#'
#' @param geneids geneids
#' @param mode go or kegg or wp
#'
#' @return enrich_res
#' @export
#'
#' @examples
#' genes=c("ASGR2","BEST1","SIGLEC16","ECRP","C1QC","TCN2","RNASE2","DYSF","C1QB","FAM20A","FCGR1A","CR1","HP","VSIG4","EGR1")
#' gene2id(genes)->geneid
#' enrich(geneid$entrezID)->GO
#' plot.enrich_res(GO$GO_res[1:10,])
#' enrichplot::cnetplot(GO$GO)
#' enrichplot::dotplot(GO$GO)
#' enrichplot::heatplot(GO$GO)
#' enrichplot::pairwise_termsim(GO$GO)%>%enrichplot::emapplot()
#' enrichplot::upsetplot(GO$GO)
enrich<-function(geneids,mode="go"){
  lib_ps("clusterProfiler","org.Hs.eg.db")
  if(mode=="go"){
    GO=enrichGO(gene = geneids,
              OrgDb = org.Hs.eg.db, # 参考基因组
              pvalueCutoff =0.05,	# P值阈值
              qvalueCutoff = 0.05,	# qvalue是P值的校正值
              ont="all",	# 主要的分为三种，三个层面来阐述基因功能，生物学过程（BP），细胞组分（CC），分子功能（MF）
              readable =T)	# 是否将基因ID转换为基因名
  }
  if(mode=="kegg")GO<-enrichKEGG(gene = geneids,keyType = "kegg",organism= "human", qvalueCutoff = 1, pvalueCutoff=1)
  if(mode=="wp") GO<-enrichWP(geneids,"Homo sapiens", qvalueCutoff = 1, pvalueCutoff=1)

  # 强制转换为数据框
  GO_all=as.data.frame(GO)
  # 筛选显著富集的数据
  GO_res<-GO_all[(GO_all$qvalue<0.05 & GO_all$p.adjust<0.05),]
  # 保存数据
  # write.table(GO,file="GOterm.csv",sep=",",quote=F,row.names = F)
  class(GO_res)<-c("enrich_res",class(GO_res))
  return(list(GO=GO,GO_res=GO_res))
}

#' Plot
#'
#' @param GO enrich_res object
#' @param mode mode
#' @param str_width default: 50
#'
#' @return ggplot
#' @exportS3Method
#'
plot.enrich_res<-function(GO,mode=1,str_width=50){
  #经典图
  if(mode==1){p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Count, fill=p.adjust))+
    geom_bar(stat = "identity",width=0.7)+####柱子宽度
    scale_fill_gradient(low = "red",high ="blue",limits=c(0,0.05))+#颜色自己可以换
    labs(x = "Gene numbers")+
    scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
    theme_bw()}

  if(mode==2){p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Count, fill=ONTOLOGY))+
    geom_bar(stat = "identity",width=0.7)+####柱子宽度
    scale_fill_manual(values = c("#66C3A5", "#8DA1CB", "#FD8D62")) + ###颜色
    labs(x = "Gene numbers")+
    scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
    theme_bw()}

  if(length(grep("^GO",GO$ID))>0)p=p+labs(title = "GO Enrichment",y = "GO Term")
  if(length(grep("^hsa",GO$ID))>0)p=p+labs(title = "KEGG Pathways Enrichment",y = "Pathway")
  if(length(grep("^WP",GO$ID))>0)p=p+labs(title = "WiKi Pathways Enrichment",y = "Pathway")
  p
}

#GSEA
if(F){
  #differ_res%>%filter(padj<0.05,abs(log2FoldChange)>1)%>%arrange(-log2FoldChange)->genes
  #rownames(genes)<-genes$tax
  #gene2id(genes$tax)->rt
  #genelist_sort=genes[rt$gene,"log2FoldChange"]
  #names(genelist_sort)=rt$entrezID

  genes=c("ASGR2","BEST1","SIGLEC16","ECRP","C1QC","TCN2","RNASE2","DYSF","C1QB","FAM20A","FCGR1A","CR1","HP","VSIG4","EGR1")
  gene2id(genes)->geneid
  genelist_sort=rnorm(14,5,3)#log2FoldChange
  names(genelist_sort)<-geneid$entrezID
  sort(genelist_sort,decreasing = T)->genelist_sort
  #https://www.jianshu.com/p/00792ef60c0d
  lib_ps("clusterProfiler","org.Hs.eg.db")
  go <- gseGO(genelist_sort,
                 ont = "ALL",
                 OrgDb = org.Hs.eg.db,
                 minGSSize    = 10,  #设置基因集范围
                 maxGSSize = 500,
                 pvalueCutoff = 1)

  go.df=as.data.frame(go)

  gseaplot(go,geneSetID = go.df$ID[1])
  enrichplot::ridgeplot(go,10,label_format =30)
}


