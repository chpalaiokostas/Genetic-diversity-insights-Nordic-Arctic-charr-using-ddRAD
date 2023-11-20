library(adegenet)
library(tidyverse)
library(reshape2)
library(vcfR)

setwd("~/Desktop/Arctic_charr_Scandinavia_genetic_diversity/")

ac_vcf <- read.vcfR("ac_scandinavia_filt_thin.vcf")

pop_info <- read.table(file="ac_scandinavia_map_filt.txt",
                       header=T,sep="\t",stringsAsFactors = F)

pop_info <- pop_info[pop_info$id %in% colnames(ac_vcf@gt),]

ac <- vcfR2genlight(ac_vcf)
pop(ac) <- pop_info$pop
ac

# allele freq 
myFreq <- glMean(ac)
myFreq <- ifelse(myFreq > 0.5, myFreq - 0.5, myFreq)
#myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)

#PCA species
pca <- glPca(ac)
myCol <- colorplot(pca$scores,pca$scores, transp=TRUE, cex=1.5)
abline(h=0,v=0, col="grey")

# estimate explained variance for first two PC
pca$eig[1]/sum(pca$eig)
pca$eig[2]/sum(pca$eig)

PC <- data.frame(pca$scores)

p1 <- ggplot(data=PC,mapping=aes(x=PC[,1],y=PC[,2])) +
  geom_point(aes(col=ac$pop),alpha=0.8,size=3) +
  labs(x="PC1 (31% explained variance)",y="PC2 (15% explained variance)",color="population") +
  #scale_colour_viridis_d() +
  scale_color_manual(values = c("Finland" = "deepskyblue3",
                                "Norway"="indianred3",
                                "Sweden"="yellow2")) +
  theme(plot.title=element_text(size=18,hjust=0.5),
        axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16))

p1 + theme_classic()

#optimal param for DAPC 

grp <- find.clusters(ac, choose.n.clust = F, criterion = "diffNgroup")

grp$size
ind_groups <- as.data.frame.matrix(table(pop(ac), grp$grp))
pop_groups <- as.data.frame.matrix(table(pop_info[,2], grp$grp))
ind_groups <- rownames_to_column(ind_groups)
colnames(ind_groups)[1] <- "Population"
ind_groups

dapc_pop <- dapc(ac,grp$grp)

## structure like plot
dapc.results <- as.data.frame(dapc_pop$posterior)
dapc.results$pop <- pop(ac)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_cluster","Posterior_membership_probability")

dapc.results <- dapc.results %>% arrange(Original_Pop)

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_cluster)) +
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = c("1" = "indianred3",
                               "2"="deepskyblue3",
                               "3"="yellow2")) +
  #scale_fill_manual(values = cols) +
  facet_grid(~Original_Pop, scales = "free") +
  labs(x="",y="Posterior membership probability",fill="Assigned\ncluster") +
  theme(axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=14),
        axis.ticks.x=element_blank()) 
p

# population prediction using CV
# 70% training - 30% validation
set.seed(19062020)
cv_results <- vector(length=100)
cv_results_inv <- data.frame()
for(i in 1:100) {
  kept.id <- unlist(tapply(1:nInd(ac), pop(ac),
                           function(x) sample_frac(as_tibble(x), size=0.7,replace=FALSE)))
  x.sup <- ac[-kept.id]
  x <- ac[kept.id]
  nInd(x)
  dapc2 <- dapc(x,n.pca=2,n.da=2)
  pred.sup <- predict.dapc(dapc2, newdata=x.sup)
  cv_results_inv <- rbind(cv_results_inv,cbind(as.character(pred.sup$assign),as.character(pop(x.sup))))
  table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))
  round(pred.sup$posterior,3)
  cv_results[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
}

colnames(cv_results_inv) <- c("Predictions","Labels")

write.table(cv_results_inv,file = "ac_scandinavia_cv_results.txt",quote = F, row.names = F, sep = "\t")

