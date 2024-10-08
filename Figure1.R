library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)

library(viridis)


# read in bulk RNA-seq data ------------------------------------------------

PBRAclusters<-read.delim("input_data/PBandRA_clusters.tsv")
RAvControl<-read.delim('input_data/RAvControl_DESeq2output.csv',sep=',')
PBvControl<-read.delim('input_data/PBvControl_DESeq2output.csv',sep=',')
PBRAvControl<-read.delim('input_data/PBandRAvControl_DESeq2output.csv',sep=',')

# box plots ---------------------------------------------------------------
RAvControl$group<-rep('RA v DMSO',nrow(RAvControl))
PBvControl$group<-rep('PB v DMSO',nrow(PBvControl))
PBRAvControl$group<-rep('PB+RA v DMSO',nrow(PBRAvControl))

alldata<-rbind(RAvControl,PBvControl,PBRAvControl)
alldata$group<-factor(alldata$group,levels=c('RA v DMSO','PB v DMSO','PB+RA v DMSO'))
alldata<-merge(alldata,PBRAclusters,by='gene')
alldata$cluster<-factor(alldata$cluster)

f1<-ggplot(alldata,aes(x=cluster,y=log2FoldChange,fill=group))+geom_boxplot(outlier.size=0.2)+
  theme_classic()+scale_fill_manual(values=c('royalblue3','red3','purple'))+coord_trans(ylim=c(-5,5))+
  guides(fill=guide_legend(title="RNA-seq comparison")) + xlab('cluster') +
  theme(axis.title = element_text(size=14),axis.text=element_text(size=12))
ggsave('output/PBRA_boxplots.png',f1,width=6,height=4)


# adrenal medulla module scores --------------------------------------------
seurat_medulla<-readRDS('scRNA_seq/adrenal_medulla_Seurat.RDS')
seurat_medulla <- UpdateSeuratObject(object = seurat_medulla)

dim(seurat_medulla)
dp<-DimPlot(seurat_medulla, reduction = "umap")
ggsave('output/adrenal_medulla_umap.png',dp,width=7.5,height=5.1)

for(cl in unique(PBRAclusters$cluster)){
  print(cl)
  signature_gene_list <-list(subset(PBRAclusters,cluster==cl)$gene)
  print(length(signature_gene_list[[1]]))
  
  scores <- AddModuleScore(object =seurat_medulla, features = signature_gene_list, name = cl)
  plot<-FeaturePlot(object = scores, features = paste0(cl,"1"), order=T)+scale_color_viridis(discrete = FALSE, option="turbo")
  plot<-plot+ggtitle(paste0('cluster ',cl))
  ggsave(filename=paste0('output/adrenal_medulla_cluster',cl,'.png'),plot=plot,width=5.5,height=5)
}


