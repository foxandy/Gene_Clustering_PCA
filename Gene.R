gene<-read.csv("gene.csv")
###data preparation
#standardize the data
Zgene<-scale(gene[,1:1000])
#already scaled
fit = prcomp(Zgene,scale=F)
#scree plot
plot(fit)
summary(fit)
#check with the covariance matrix
geneCov<-cov(Zgene[,1:1000])
fit2=eigen(geneCov)
(fit2$values[1]+fit2$values[2])/(sum(fit2$values))
###PCA
sick<-gene$sick
genePC<-as.data.frame(cbind(fit$x,sick))
plot(genePC$PC1,genePC$PC2,col=genePC$sick+3,main="Gene Principal Component Analysis",xlab="PC1",ylab="PC2")
###Hierarchical Clustering
geneDist<-dist(Zgene,method="euclidean")
hier=hclust(geneDist,method="complete")
plot(hier,labels=sick,main="Hierarchical Clustering of Gene Measurements")
corrDist<-as.dist(1-cor(t(gene[,1:1000])))
hierCorr<-hclust(corrDist,method="complete")
plot(hierCorr,labels=sick,main="Hierarchical Clustering using Correlation")
###K-means Clustering
fit.Kmean<-kmeans(Zgene,2,nstart=100)
summary(fit.Kmean)
plot(fit.Kmean)
K.clusters<-table(sick,fit.Kmean$cluster)
K.clusters
###Gaussian Mixtures
library(mclust)
fit.Gaussian<-Mclust(genePC[,1:2],G=2,modelNames=c("EEI"))
plot(fit.Gaussian)
fit.Gaussian$parameter$pro
fit.Gaussian$parameter$mean
fit.Gaussian$parameter$variance
table(sick,fit.Gaussian$classification)
min(genePC[1:20,1])
max(genePC[1:20,1])
min(genePC[21:40,1])
max(genePC[21:40,1])
plot(fit,main="PCA Scree Plot",xlab="PC")
