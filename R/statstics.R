#======some statstics models

#' Multiple regression/ variance decomposition analysis
#'
#' @param data dataframe
#' @param TopN give top variable importance
#' @param formula formula
#'
#' @examples
#' data(otutab)
#' cbind(a_diversity(otutab,method="shannon"),metadata[,c(2,12:15)])->a_res
#' multireg(formula = Shannon~Group*.,data = a_res)
multireg<-function(formula,data,TopN=3){
  model.frame(formula,data=data)->metatbl
  colnames(metatbl)[1:2]<-c("test_v","xGroup")
  metatbl$xGroup<-factor(metatbl$xGroup)

  lib_ps("relaimpo","aplot")
  #multi-lm and correlation
  n_env_cor=list()
  n_env_lm=list()
  for (i in levels(metatbl$xGroup)){
    print(i)
    metadata=metatbl%>%filter(xGroup==i)
    n_env_cor[[i]]=cor(metadata[,"test_v"],metadata[,-1:-2])
    mlm<-lm(test_v~.,cbind(metadata[,"test_v",drop=F],metadata[,-1:-2]))
    relaimpo::calc.relimp(mlm)->mlm_r
    n_env_lm[[i]]=mlm_r@lmg
  }

  do.call(rbind,n_env_cor)%>%as.data.frame%>%mutate(xGroup=names(n_env_cor))->n_env_cor
  do.call(rbind,n_env_lm)%>%as.data.frame%>%mutate(xGroup=names(n_env_lm))->n_env_lm

  melt(n_env_cor,id.vars = "xGroup")->n_env_cor
  melt(n_env_lm,id.vars = "xGroup")->n_env_lm
  n_env_lm%>%group_by(xGroup)%>%summarise(explained=sum(value))->n_explained
  n_env_lm%>%group_by(xGroup)%>%top_n(TopN)->n_env_lm1

  p1=ggplot(data =n_env_cor,aes(y=variable,x=xGroup))+
    geom_tile(aes(fill=value))+
    scale_fill_gradient2(name="Correlation",low = "#6D9EC1", mid = "white", high = "#E46726",midpoint = 0)+
    geom_point(data = n_env_lm1,aes(size=value),shape=21,fill=NA)+
    guides(size=guide_legend("Importance"))+labs(x=NULL,y=NULL)+theme_pubr(base_size = 14,legend = "right")

  p2=ggbarplot(n_explained,x = "xGroup",y="explained",fill = "#4EA9E6")+
    scale_y_continuous(expand = c(0,0))+labs(x=NULL,y=NULL,subtitle ="Explained variation")+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),plot.margin = unit(c(0,0,3,0),"lines"))

  p1%>%insert_top(p2,height = 0.3)
}
