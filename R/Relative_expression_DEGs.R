library(dplyr)
library(factoextra)
library(ComplexHeatmap)
library(clustree)
library(ggplot2)
temp<-list.files(path = "tsv_files/", pattern = "*_DE.tsv")
myfiles<-lapply(paste0("tsv_files/",temp),function(i){read.table(i,header = T,sep="\t")})
targets<-unlist(read.csv("csv_files/Candidate_genes.csv",header = F,sep = ","))

# Subset the DEGs from each *.tsv files and merge them
df<-list()
for (i in 1:length(temp)) {
  df[[i]] <- myfiles[[i]]
  df[[i]] <- filter(df[[i]], pval<10e-6 & abs(log2FoldChange)>2)
  df[[i]]<-select(df[[i]],gene_id,log2FoldChange)
  names(df[[i]])<-c("features",targets[i])
}
feat<-df[[1]]$features
for (i in 2:length(temp)) {
  feat<-unique(c(feat,df[[i]]$features))
}
df_filter<-list()
for (i in 1:length(temp)) {
  df_filter[[i]]<-filter(myfiles[[i]],gene_id %in% feat & pval<10e-6)
  df_filter[[i]]<-select(df_filter[[i]],gene_id,log2FoldChange)
  names(df_filter[[i]])<-c("features",targets[i])
}
df_re_exp<-df_filter[[1]]
for (i in 2:length(temp)) {
  df_re_exp<-full_join(df_re_exp,df_filter[[i]],by="features")
}
df_re_exp[is.na(df_re_exp)]<-0
df_re_exp$NonTargeting<-0
write.csv(df_re_exp,"csv_files/relative_expression.csv",row.names = F)

# Drawing the heatmap of the relative expression of each DEGs:

df_re_exp<- read.csv("csv_files/relative_expression.csv",header = T,sep = ",")

# Transpose the dataframe, set the column names with the features and the row names with groups
df_t<-t(df_re_exp[,-1])
colnames(df_t)<-df_re_exp[,1]

# Determine and visualize the optimal number of clusters using hierarchical clustering
fviz_nbclust(df_t,FUN=hcut, method = "silhouette")

# Hierarchical clustering 
d <- dist(df_t, method = "euclidean")
res.hc <- hclust(d, method = "ward.D2" )

# Get the cluster tree:
tmp <- NULL
for (i in 1:63){
  tmp[[i]] <- data.frame(cluster=cutree(as.hclust(res.hc),i))
}
df <- data.frame(tmp)
colnames(df) <- seq(1:63)
colnames(df) <- paste0("clusters",colnames(df))
df<-mutate(df,feature=rownames(df_t))

# save the cluster tree in pdf or svg:
pdf("clustree.pdf",height=50,width = 20)
clustree(df, prefix = "clusters")
dev.off()
svglite::svglite("clustree.svg",height=50,width = 20)
clustree(df, prefix = "clusters")
dev.off()


