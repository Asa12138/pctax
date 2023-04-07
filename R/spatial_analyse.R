if(F){#Spatial variables construct
library(adespatial)
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")
library(spacemakeR)
data("otutab")
#some missing values in original coordinates data, run imputation first from the RDA section
rda.data=decostand(t(otutab),'hellinger')%>%as.data.frame()
expo.xy.name=metadata[,c(9,10,1)]

colnames(expo.xy.name) = c("latitude", "longitude", "name")
expo.xy=expo.xy.name[,1:2]

geosphere::distm(expo.xy.name[,2:1])->expo.geo.xy
rownames(expo.xy.name)->rownames(expo.geo.xy)->colnames(expo.geo.xy)
expo.geo.xy=round(expo.geo.xy/1000)#distance in km


expo.vario = variogmultiv(Y = rda.data, xy = expo.xy, nclass =20)
plot(expo.vario$d, expo.vario$var, ty='b', pch=20, xlab="Distance", ylab="C(distance)") #210 is the maxmium
#first peak is
expo.pcnm = vegan::pcnm(dist(expo.xy))
min.d = expo.pcnm$threshold #minimum distance
thresh10 = seq(give.thresh(dist(expo.xy)), give.thresh(dist(expo.xy))*4, le=10)
list10nb = lapply(thresh10, dnearneigh, x=as.matrix(expo.xy), d1=0)
f2 = function(D, dmax, y) {1- (D/dmax)^y}
expo.thresh.f2 = lapply(list10nb, function(x) test.W(x, Y=rda.data, f=f2, y=2:10, dmax=max(unlist(nbdists(x, as.matrix(expo.xy)))), xy=as.matrix(expo.xy)))

expo.f2.minAIC = sapply(expo.thresh.f2, function(x) min(x$best$AICc, na.rm=TRUE))
min(expo.f2.minAIC)
nb.bestmod = which.min(expo.f2.minAIC) #6 variables
dmax.best = expo.thresh.f2[nb.bestmod][[1]]$all[1,2]
expo.MEM.champ = unlist(expo.thresh.f2[which.min(expo.f2.minAIC)], recursive = FALSE)
summary(expo.MEM.champ)
expo.MEM.champ$best$values
expo.MEM.champ$best$ord
MEMid = expo.MEM.champ$best$ord[1:which.min(expo.MEM.champ$best$AICc)]
sort(MEMid)
expo.MEM.champ$best$R2
MEM.select = expo.MEM.champ$best$vectors[, sort(c(MEMid)),drop=F]
colnames(MEM.select) = sort(MEMid)
#best R2
R2.MEMbest = expo.MEM.champ$best$R2[which.min(expo.MEM.champ$best$AICc)]
#adjusted best R2
RsquareAdj(R2.MEMbest, nrow(rda.data), length(MEMid))
par(mfrow=c(2,4))
for(i in 1:ncol(MEM.select)){ #visualizing selected MEM factors
  s.value(expo.xy, MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}

expo.MEM.rda = rda(rda.data ~., as.data.frame(MEM.select))
expo.MEM.R2a = RsquareAdj(expo.MEM.rda)$adj.r.squared
anova(expo.MEM.rda)
axes.MEM.test = anova(expo.MEM.rda, by="axis")
nb.ax = length(which(axes.MEM.test[,5] <= 0.05))

#plot maps of the significant canonocial axes
expo.MEM.axes = scores(expo.MEM.rda, choices=c(1,2), display="lc", scaling=1)
par(mfrow = c(1,2))
s.value(expo.xy, expo.MEM.axes[,1])
s.value(expo.xy, expo.MEM.axes[,2])

#visualizing top two axes (tested by axes.MEM.test)
for (i in 1:2){
  print(ggplot(expo.xy, aes(y=latitude, x=longitude, size = expo.MEM.axes[,i])) + geom_point(shape=21, alpha=0.6, fill="steel blue") + labs(size = paste0("RDA.MEM.axis ", i)))
}

#visualizing MEM.select using ggplot
for (i in 1:ncol(MEM.select)){
  print(ggplot(expo.xy, aes(y=latitude, x=longitude, size = MEM.select[,i])) + geom_point(shape=21, alpha=0.6, fill="steel blue") + labs(size = paste0("MEM variable ", i)))
}

#testing location~season randomness
library(nnet)
locaseason = cbind(as.data.frame(plot.MEM.select), season = scores_single_complete$season)
test <- multinom(season ~ ., data = locaseason)
summary(test)
}
