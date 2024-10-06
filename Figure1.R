library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)
#library(tidyr)
library(viridis)


# clusters data -----------------------------------------------------------

PBRAclusters<-read.delim('scRNA_seq/output_data/PBandRA_cluster_k5_complete.txt',sep=' ')
group_names<-c('PB+RA upregulated cluster','PB+RA and PB regulated','all downregulated','RA upregulated','PB upregulated')
relabel_clusters<-c(4,5,1,2,3)
for(i in 1:nrow(PBRAclusters)){
  PBRAclusters[i,'condition']<-group_names[PBRAclusters[i,'cluster']]
  PBRAclusters[i,'cluster_v2']<-relabel_clusters[PBRAclusters[i,'cluster']]
}


# box plots ---------------------------------------------------------------

RAvControl<-read.delim('scRNA_seq/output_data/RAvControl_DESeq2output.csv',sep=',')
PBvControl<-read.delim('scRNA_seq/output_data/PBvControl_DESeq2output.csv',sep=',')
PBRAvControl<-read.delim('scRNA_seq/output_data/PBandRAvControl_DESeq2output.csv',sep=',')

RAvControl$group<-rep('RA v DMSO',nrow(RAvControl))
PBvControl$group<-rep('PB v DMSO',nrow(PBvControl))
PBRAvControl$group<-rep('PB+RA v DMSO',nrow(PBRAvControl))

alldata<-rbind(RAvControl,PBvControl,PBRAvControl)
alldata$group<-factor(alldata$group,levels=c('RA v DMSO','PB v DMSO','PB+RA v DMSO'))
alldata<-merge(alldata,PBRAclusters,by='gene')
alldata$cluster_v2<-factor(alldata$cluster_v2)

f1<-ggplot(alldata,aes(x=cluster_v2,y=log2FoldChange,fill=group))+geom_boxplot(outlier.size=0.2)+
  theme_classic()+scale_fill_manual(values=c('royalblue3','red3','purple'))+coord_trans(ylim=c(-5,5))+
  guides(fill=guide_legend(title="RNA-seq comparison")) + xlab('cluster') +
  theme(axis.title = element_text(size=14),axis.text=element_text(size=12))
ggsave('new_plots/PBRA_boxplots.png',f1,width=6,height=4)


# adrenal medulla ---------------------------------------------------------
seurat_medulla<-readRDS('scRNA_seq/adrenal_medulla_Seurat.RDS')
seurat_medulla <- UpdateSeuratObject(object = seurat_medulla)

dim(seurat_medulla)
dp<-DimPlot(seurat_medulla, reduction = "umap")
ggsave('new_plots/adrenal_medulla_umap.png',dp,width=7.5,height=5.1)
dp<-DimPlot(seurat_medulla, reduction = "umap",label=TRUE,label.size=2)+
  theme(legend.position = "none")
ggsave('new_plots/adrenal_medulla_umap_v2.png',dp,width=5,height=5)
fp<-FeaturePlot(object = seurat_medulla, features = 'nFeature_RNA')
ggsave('new_plots/adrenal_medulla_nFeature_RNA.png',fp,width=5,height=5)
fp<-FeaturePlot(object = seurat_medulla, features = 'nCount_RNA')
ggsave('new_plots/adrenal_medulla_nCount_RNA.png',fp,width=5,height=5)


meta <- seurat_medulla@meta.data
dim(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

for(cl in unique(PBRAclusters$cluster_v2)){
  print(cl)
  signature_gene_list <-list(subset(PBRAclusters,cluster_v2==cl)$gene)
  print(length(signature_gene_list[[1]]))
  
  scores <- AddModuleScore(object =seurat_medulla, features = signature_gene_list, name = cl)
  plot<-FeaturePlot(object = scores, features = paste0(cl,"1"), order=T)+scale_color_viridis(discrete = FALSE, option="turbo")
  plot<-plot+ggtitle(paste0('cluster ',cl))
  ggsave(filename=paste0('new_plots/adrenal_medulla_cluster',cl,'.png'),plot=plot,width=5.5,height=5)
}

# each plot separately with each cluster and sorted naming based on grouping



adrenal_medulla_genes<-Features(seurat_medulla)
length(adrenal_medulla_genes)

for(cl in unique(PBRAclusters$cluster)){
  print(cl)
  genes<-subset(PBRAclusters,cluster==cl)$gene
  print(length(genes))
  print(length(genes[genes %in% adrenal_medulla_genes]))
}
