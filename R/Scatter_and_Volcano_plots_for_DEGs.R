library(EnhancedVolcano)
library(dplyr)
library(magrittr)
library(grid)
library(gridExtra)
library(ggplot2)

temp<-list.files(path = "tsv_files/", pattern = "*_DE.tsv")
myfiles<-lapply(paste0("tsv_files/",temp),function(i){read.table(i,header = T,sep="\t")})
targets<-unlist(read.csv("csv_files/Candidate_genes.csv",header = F,sep = ","))

# function of drawing plots
dual_plots<-function(df,target) {
  df<-mutate(df,color=case_when(pval<10e-6 & log2FoldChange>2 ~ "Upregulated",
                                pval<10e-6 & log2FoldChange<(-2) ~ "Downregulated",
                                TRUE ~ "Black"))
  n_up<-nrow(filter(df,color=="Upregulated"))
  n_down<-nrow(filter(df,color=="Downregulated"))
  
  m<-ggplot(df,aes(y=log2(test_mean+1),x=log2(control_mean+1)))+
    geom_point(aes(color=as.factor(color)),alpha=0.44)+
    scale_color_manual(breaks=c("Upregulated","Downregulated","Black"),
                       limits=c("Upregulated","Downregulated","Black"),
                       values = c("green","red","grey40"))+
    geom_text_repel(data=df %>% filter(pval<10e-6 & abs(log2FoldChange)>2),
                    aes(label=gene_id),hjust=0,vjust=0)+theme_bw()+
    theme(panel.grid.major = element_line(size=1.5,colour = "grey96"),
          panel.grid.minor = element_line(size=0.2,colour = "grey96"),
          panel.border = element_blank(),
          axis.line = element_line(size=1,colour = "black"),
          legend.position = "none",
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size=12))+
    ylab(paste0("Log2(",target,"_Knockedout cells+1)",sep = ""))+
    xlab("Log2(Nontargeting_Control cells+1)")
  
  n<-ggplot(df,aes(y=-log10(pval),x=log2FoldChange))+
    geom_point(aes(color=as.factor(color)),alpha=0.44)+
    xlim(-5,5)+
    geom_hline(yintercept = -log10(10e-6),linetype="dashed")+
    geom_vline(xintercept = -2,linetype="dashed")+
    geom_vline(xintercept = 2,linetype="dashed")+
    geom_text_repel(data=df %>% filter(pval<10e-6 & abs(log2FoldChange)>2),
                    aes(label=gene_id),hjust=0,vjust=0)+theme_bw()+
    theme(panel.grid.major = element_line(size=1.5,colour = "grey96"),
          panel.grid.minor = element_line(size=0.2,colour = "grey96"),
          panel.border = element_blank(),
          axis.line = element_line(size=1,colour = "black"),
          legend.position = "right",
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size=12))+
    ylab("-Log10(p value)")+
    xlab("Log2 Fold Change")+
    scale_color_manual(name=NULL,
                       breaks=c("Upregulated","Downregulated","Black"),
                       limits=c("Upregulated","Downregulated","Black"),
                       values=c("green","red","grey40"),
                       labels=c("Upregulated DEGs","Downregulated DEGs","NS"))
  
  print(grid.arrange(m,n,nrow=1,widths=4:5, 
                  top=textGrob(paste(target,"vs Nontargeting Controls"),
                               x=0.45,
                               gp=gpar(fontsize=14,fontfact="bold"),
                               vjust = 0.6),
                  bottom=textGrob(paste("Up_DEGs: ",n_up,", Down_DEGs: ",n_down,
                                        ", pvalue: 10e-6, Log2_FC: 2", sep = ""),
                                  x=0.85,
                                  gp=gpar(fontsize=12))))
  
}

pdf("dual_plots.pdf",height = 7.5,width=15,bg="transparent")
for (i in 1:length(temp)) {dual_plots(myfiles[[i]],targets[i])} 
dev.off()

