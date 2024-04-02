library(StatMatch) 
library(cluster)
library(factoextra)
library(clusterSim)
library(vegan)
library(biotools)
library(ggplot2)
library(dplyr)
library(MVN)

eco <- read.csv("C:\\Users\\Sam\\Documents\\Spring 24\\multivariate\\task 2\\ecosystems.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE ,dec=",")



#scale the data from 0-1
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
eco <- as.data.frame(lapply(eco, min_max_scale))

#use mahalanobis distance
#DM<-as.dist(mahalanobis.dist(eco)) 
#d.eco = as.matrix(DM)

#use hellinger distance
#data.hel = decostand(eco, method="hel")
#d.eco = dist(data.hel)


#use BC distance
data.hel = vegdist(eco, method="bray")
d.eco = dist(as.matrix(data.hel))



cluster.H.AV <- hclust(DM, method="ward.D2")
plot(cluster.H.AV,hang = -1, cex=0.7)




PseudoF<-data.frame()
for (centers in c(2,3,4,5,6,7,8,9,10,11)){
  km.hel<-pam(d.eco,centers,diss=T) # Partitioning Around Medoids
  PseudoF.hel<-index.G1(x=eco, cl=km.hel$cluster, d=d.eco, centrotypes = "medoids")
  PseudoF[(centers-1),1]<-centers
  PseudoF[(centers-1),2]<-PseudoF.hel
  5
}
PseudoF
plot(PseudoF[,1],PseudoF[,2],type="b",xlab="k",ylab="PseudoF")




sil<-fviz_nbclust(eco, FUNcluster = pam, diss=NULL, method = "sil")
sil



k <- 2
pam.k <- pam(d.eco,k,diss=T)
plot(pam.k)

mds.pam <- cmdscale(d.eco, k= k, eig=TRUE)
mds.pam$GOF

clusplot(mds.pam$points, pam.k$clustering, lines=0, shade=T, main="PAM", xlab="Axis 1", ylab="Axis 2",color=T)

boxM(eco,pam.k$clustering)



#display boxplot of group columns
pam.k <- pam(d.eco, k, diss=T)
eco$cluster <- pam.k$clustering
plots <- list()
for (i in 1:(ncol(eco) - 1)) {
  column_name <- names(eco)[i]
  p <- ggplot(eco, aes_string(x = "factor(cluster)", y = column_name, fill = "factor(cluster)")) +
    geom_boxplot() +
    labs(title = paste(column_name), 
         x = "Cluster", 
         y = column_name)
  plots[[i]] <- p
}
library(gridExtra)
do.call(grid.arrange, c(plots, ncol=3))


g1<-eco[pam.k$clustering==1, -which(names(eco) == "cluster")]
dim(g1)
det(cov(g1))
mvn(g1,mvnTest="mardia",multivariatePlot="qq")


g2<-eco[pam.k$clustering==2, -which(names(eco) == "cluster")]
dim(g2)
det(cov(g2))
mvn(g2,mvnTest="mardia",multivariatePlot="qq")



boxM(eco[, -which(names(eco) == "cluster")],pam.k$clustering)

eco.manova<-manova(cbind(Colif_total, Colif_fecal, Estrep_fecal, Cont_mineral, Conductivitat, Solids_susp, DQO_M) ~ pam.k$clustering, 
                      data=eco[, -which(names(eco) == "cluster")])
summary(eco.manova, test = "Pillai")







new_data <- data.frame(
  Colif_total = c(202, 204),  
  Colif_fecal = c(193, 198),  
  Estrep_fecal = c(97, 100),  
  Cont_mineral = c(42, 37),  
  Conductivitat = c(21, 21),  
  Solids_susp = c(1026, 979),  
  DQO_M = c(714, 714)   
)

eco2 <- rbind(eco[, -which(names(eco) == "cluster")], new_data)
eco2 <- as.data.frame(lapply(eco2, min_max_scale))
new_dat <- eco2[(nrow(eco2) - 1):nrow(eco2), ]


d11 <- vegdist(rbind(new_dat[1,], eco[pam.k$medoids[1], c(1:7)]), method = "bray")
d12 <-vegdist(rbind(new_dat[1,], eco[pam.k$medoids[2], c(1:7)]), method = "bray")
d21 <- vegdist(rbind(new_dat[2,], eco[pam.k$medoids[1], c(1:7)]), method = "bray")
d22 <- vegdist(rbind(new_dat[2,], eco[pam.k$medoids[2], c(1:7)]), method = "bray")

distance_matrix <-  matrix(c(d11, d12, d21, d22), nrow = 2, byrow = TRUE)
print(distance_matrix)
