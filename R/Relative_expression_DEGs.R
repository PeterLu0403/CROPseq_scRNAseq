library(dplyr)
library(factoextra)
library(ComplexHeatmap)
library(clustree)
library(ggplot2)
library(dendextend)
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
write.csv(df_re_exp,"outs/relative_expression.csv",row.names = F)

# Do clustering:
df_re_exp<- read.csv("outs/relative_expression.csv",header = T,sep = ",")
# Transpose the data frame, set the column names with the features and the row names with groups
df_t<-t(df_re_exp[,-1])
colnames(df_t)<-df_re_exp[,1]
# Determine and visualize the optimal number of clusters using hierarchical clustering
fviz_nbclust(df_t,FUN=hcut, method = "silhouette")
# Hierarchical clustering 
d <- dist(df_t, method = "euclidean")
res.hc <- hclust(d, method = "ward.D2" )
# Drawing a heatmap of the relative expression
cut_row<-data.frame(cluster=cutree(as.hclust(res.hc),2))
cut_row<-mutate(cut_row,feature=rownames(cut_row))
pdf("outs/relative_expression.pdf",height=15,width=25,bg="transparent")
Heatmap(df_t, show_row_names = T,show_column_names = F, 
        row_split = as.factor(cut_row$cluster),
        cluster_columns = T,
        heatmap_legend_param = list(title = "Log2FC"))
dev.off()
# Get the cluster tree:
tmp <- NULL
for (i in 1:63){
  tmp[[i]] <- data.frame(cluster=cutree(as.hclust(res.hc),i))
}
df_tree<- data.frame(tmp)
colnames(df_tree) <- seq(1:63)
colnames(df_tree) <- paste0("clusters",colnames(df_tree))
df_tree<-mutate(df_tree,feature=rownames(df_t))
# save the tree file as a csv file:
write.csv(df_tree,"outs/clustree.csv",row.names = F)
# save the cluster tree in pdf or svg:
pdf("outs/clustree.pdf",height=50,width = 20)
clustree(df_tree, prefix = "clusters",edge_arrow_ends = "first")
dev.off()
svglite::svglite("outs/clustree.svg",height=50,width = 20)
clustree(df_tree, prefix = "clusters",edge_arrow_ends = "first")
dev.off()



