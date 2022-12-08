#======some statstics models

#' Multiple regression/ variance decomposition analysis
#'
#'
#'
#'
if(F){
lib_ps("relaimpo")

data(otutab)
a_diversity(otutab)->a_res
mlm<-lm(Shannon~.,cbind(a_res[,"Shannon",drop=F],metadata[,12:19]))
summary(mlm)

calc.relimp(mlm)->mlm_r

data.frame(factor=names(mlm_r@lmg),`relative importance`=mlm_r@lmg,check.names = F)%>%
  ggdotchart(.,"factor","relative importance",sorting = "ascending",
             dot.size = 7, color="#88ABEB",
             add = "segment",add.params = list(color = "lightgray", size = 1.5))
}
