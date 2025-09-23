#!/home/chengennong/tools/mambaforge/envs/cnmf/bin/R

library(factoextra)
library(dplyr)
library(FactoMineR)
meta_df <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% select(Sample,Group,Donor,Cell_collection_timepoint)
meta_df

pca_assay <- readRDS("/home/chengennong/ylp/aav/rnaseq/plot/vst_pca_assay.rds")
pca.res <- PCA(pca_assay, graph = F, scale.unit = T) # 简简单单1行代码实现主成分分析
pca.res

pca_assay[1:5,1:5]
pca_assay$Sample <- rownames(pca_assay)
pca_assay <- merge(pca_assay, meta_df,by="Sample")
dim(pca_assay)
#p <- fviz_pca_ind(res.pca, label="none", habillage=meta_df$Condition,
#             addEllipses=TRUE, ellipse.level=0.95)+
colnames(pca.res$eig)
#p <- fviz_pca_ind(pca.res,
#             col.ind = meta_df$Cell_collection_timepoint,    # 用 condition 着色
#             #col.ind = meta_df$Donor,    # 用 condition 着色
#             #shape.ind = meta_df$Group,      # 用 time 设置形状
#             palette = "Set1",             
#             #addEllipses = TRUE,           
#             legend.title = "Condition and Time", 
#             #habillage=pca_assay$Group,
#             label="none",
#             pointsize = 1)
#
#ggsave("/home/chengennong/ylp/aav/rnaseq/plot/pca_vst.pdf", width = 6, height = 4)

pca_df <- data.frame(pca.res$ind$coord) 
pca_var <- pca.res$eig[,2]

pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, meta_df,by="Sample")
p <- ggplot(pca_df, aes(x=Dim.1, y=Dim.2,shape=paste(Donor, Group),fill=Donor, colour=Cell_collection_timepoint)) + 
  geom_point(size=3) +
  theme_bw() +
  #scale_shape_manual(values = c(15,17,0,2),labels=c("Salt in Fanggan","Salt in Panjin","Fresh in Fanggan","Fresh in Panjin"))+
  stat_ellipse(aes(shape = NULL,linetype=Donor))
p <- ggplot(pca_df, aes(x=Dim.1, y=Dim.2,shape=paste(Donor),fill=Donor, colour=Group,alpha=Cell_collection_timepoint)) + 
  geom_point(size=3) +
  theme_bw() +
  #scale_shape_manual(values = c(15,17,0,2),labels=c("Salt in Fanggan","Salt in Panjin","Fresh in Fanggan","Fresh in Panjin"))+
  stat_ellipse(aes(color = NULL,linetype=Donor, alpha = Cell_collection_timepoint))+
  scale_alpha_manual(name="Time point",values = c(0.4,1), labels = c("16 hours after elec", "7 days after elec"))+
    scale_color_manual(name="Condition",values=c("black","#999999", "#E69F00" ), labels = c("elec Cas9", "elec Cas9 +AAV6", "elec Cas9 +AAV.APR31"))+
    labs(shape = "Donor", x= paste0("PC1 (",round(pca_var[1],2),"%)"), y= paste0("PC2 (",round(pca_var[2],2),"%)"))
ggsave("/home/chengennong/ylp/aav/rnaseq/plot/pca_vst.pdf", width = 5.8, height = 4)
