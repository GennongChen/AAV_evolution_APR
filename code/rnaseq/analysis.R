#!/usr/bin/R



#meta data ####

meta_df <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% select(Sample,Group,Donor,Cell_collection_timepoint)
rownames(meta_df) <- meta_df$Sample
meta_df
#load data ####
library(tidyverse)
count_file <- "/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt"
counts_df <- read.table(count_file, header=T,check.names=F) %>% column_to_rownames("Geneid")
colnames(counts_df)
colnames_count <- sub("/.*", "", colnames(counts_df))
if (!setequal(colnames_count, meta_df$Sample)) {
  stop("not equal to meta_df")
}
colnames(counts_df) <- colnames_count
counts_df <- counts_df[,meta_df$Sample]
head(counts_df)
#deseq2
library(DESeq2)  #加载包
countData <- counts_df %>% as.matrix
#colData <- data.frame(row.names=meta_df$Sample, Condition=meta_df$Condition)
colData <- meta_df #%>% select(Group,Cell_collection_timepoint)

#正式构建dds矩阵
#dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ Group + Cell_collection_timepoint + Group:Cell_collection_timepoint)

countsThreshold <- 5
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= countsThreshold) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
head(dds)  #查看一下构建好的矩阵
dds <- DESeq(dds)   #对原始dds进行normalize
vst_data <- vst(dds, blind = FALSE)
pca_assay <- as.data.frame(t(assay(vst_data)))
saveRDS(pca_assay, "/home/chengennong/ylp/aav/rnaseq/plot/vst_pca_assay.rds")
## 查看结果的名称，本次实验中是 "Intercept"，"condition_akap95_vs_control"
#resultsNames(dds)
## 将结果用results()函数来获取，赋值给res变量
#res <- results(dds)
## summary一下，看一下结果的概要信息
#summary(res)
## 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
#table(res$padj<0.05) #取P值小于0.05的结果
#res <- res[order(res$padj),]
#diff_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
#diff_gene_deseq <- row.names(diff_gene_deseq2)
#resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
##得到csv格式的差异表达分析结果
#write.csv(resdata,file= "/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/Time_Tech.csv",row.names = F)
#write.csv(diff_gene_deseq2,file= "/home/chengennong/ylp/aav/rnaseq/plot/Time_Tech_DEG_LFC05.csv",row.names = T)

# sample dist ####
sampleDists <- dist(t(assay(vst_data)))
library(pheatmap)
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap::pheatmap(sampleDistMatrix,
         #show_colnames = TRUE,
         show_rownames = TRUE,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         width = 7,
         height = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot/p_sample_dist.pdf")

#heatmap ####
fea<-c("LDHA","PDCD1","HAVCR2","LAG3","TOX","NR4A2","S1PR1","TCF7","SELL","IL7R","FOXO1","GZMK","GZMB","GNLY")
counts_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt",check.names=F, header=T) 
length_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/counts.txt",check.names=F, header=T) %>% dplyr::select(1,6)
source("/home/chengennong/ylp/others/scripts/rnaseq.R")
tpm_value <- counts_to_tpm(counts_df0, length_df0)

colnames_count <- sub("/.*", "", colnames(tpm_value))
colnames(tpm_value) <- colnames_count
tpm_value <- tpm_value[,meta_df$Sample]
head(tpm_value)

cols <- colorRampPalette(colors = c("#8DD3C7", "white", "#C54B8C"))(50)
mat  <- (log(tpm_value+1))[c(fea),]
#mat  <- assay(vst_data)[c("LDHA","PDCD1","HAVCR2","LAG3","TOX","NR5A2","S1PR1","TCF7","SELL","IL7R","FOXO1"),]

mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))
dim(t(anno))
anno <- as.data.frame(colData(vst_data)[, c("Cell_collection_timepoint","Donor","Group")]) %>% arrange(Cell_collection_timepoint)
mat <- mat[rownames(anno),]
pheatmap::pheatmap(t(mat), annotation_col = anno,
         col = cols,
         width = 12,
         cluster_rows = F,
         cluster_cols = F,
         height = 7,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot/p_heatmap_gene.pdf")

