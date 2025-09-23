#!/usr/bin/R


library(tidyverse)

#meta data ####
#Donors=c("Donor 300F", "Donor 1100F")
Time_ps=c("7 Days after elec", "16hours after elec")
Vectors=c("HONG31", "AAV6")

#Donor_i=Donors[1]
#Time_p_i=Time_ps[1]
#Vector_i=Vectors[1]

test_deg_in <- function(Vector_i) {
    #load data ####
    counts_df <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt", header=T,check.names = F) %>% column_to_rownames("Geneid")
    colnames(counts_df)
    meta_df <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% dplyr::select(Sample,Group,Donor,Cell_collection_timepoint)
    counts_df <- counts_df #%>% dplyr::select(matches(paste(meta_df$Sample, collapse = "|")))
    colnames(counts_df) <- c(meta_df$Sample)
    
    rownames(meta_df) <- meta_df$Sample
    meta_df
    #names(counts_df)
    head(counts_df)
    #deseq2
    library(DESeq2)  #加载包
    countData <- counts_df %>% as.matrix
    #colData <- data.frame(row.names=meta_df$Sample, Condition=meta_df$Condition)
    colData <- meta_df #%>% select(Group,Cell_collection_timepoint)
    meta_df2 <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% dplyr::select(Sample,Group,Donor,Cell_collection_timepoint) %>% 
      filter((str_detect(Group,Vector_i)))
    #正式构建dds矩阵
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ Donor + Group)
    dds <- estimateSizeFactors(dds)
    Nor_countsThreshold <- 2
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds, normalized=TRUE) >= Nor_countsThreshold) >= smallestGroupSize
    #keep2 <- c("3F-31-7-1","11F-0-7-3")
    dds <- dds[keep,]
    counts_df2 <- counts_df[,meta_df2$Sample]
    nrow(dds)
    head(dds)  #查看一下构建好的矩阵
    dds <- DESeqDataSetFromMatrix(countData=counts_df2, colData=meta_df2, design= ~ Donor + Cell_collection_timepoint)
    dds <- DESeq(dds)   #对原始dds进行normalize
    vst_data <- vst(dds, blind = FALSE)
    pca_assay <- as.data.frame(t(assay(vst_data)))
    #saveRDS(pca_assay, "/home/chengennong/ylp/aav/rnaseq/plot_time/vst_pca_assay_tech.rds")
    # 查看结果的名称，本次实验中是 "Intercept"，"condition_akap95_vs_control"
    resultsNames(dds)
    # 将结果用results()函数来获取，赋值给res变量
    res <- results(dds)
    # summary一下，看一下结果的概要信息
    summary(res)
    # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
    table(res$padj<0.05) #取P值小于0.05的结果
    res <- res[order(res$padj),]
    resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)

    #head(resdata)
    resdata_filtered_f <- resdata %>% select(Row.names:padj) %>% 
      mutate(Vector=Vector_i) 
        return(resdata_filtered_f)
}

deg_result_raw <-data.frame()
  for (p3 in Vectors) {
    print(paste(p3,sep="_"))
    deg_result_raw <- bind_rows(deg_result_raw, test_deg_in(p3))
  }

deg_result_raw%>% head

lfc_thre=1;qval_thre=0.05
deg_result <- deg_result_raw %>% #filter(group == AL[1] & avgExpr!=0) %>% 
      mutate(Label = ifelse((log2FoldChange) > lfc_thre & padj < qval_thre, "Up", 
                            ifelse((log2FoldChange) < -lfc_thre & padj < qval_thre, "Down","None"))) %>% 
      mutate(Label=ifelse(is.na(Label), "None", Label)) %>% 
      mutate(Label=factor(Label, levels= c("Up", "Down", "None"))) %>% 
      arrange(-abs(log2FoldChange), .by_group=TRUE) %>% 
      dplyr::select(Row.names, baseMean, log2FoldChange,padj,Label,Vector) 

deg_result %>% group_by(Label,Vector) %>% summarise(count=n())

deg_result %>% filter(is.na(Label)) %>% head
deg_result%>%filter(Row.names=="WT1")%>%as.data.frame
deg_result%>%filter(Row.names=="IGHG1")%>%as.data.frame


#Cor

deg_result%>%filter(Row.names=="IGHA2")%>%as.data.frame
wide_data%>%as.data.frame%>%group_by( Label_AAV6, Label_HONG31)%>%summarise(count=n())

wide_data <- deg_result %>% #dplyr::select(-Comb) %>%
  pivot_wider(
    names_from = Vector,           # 将 Vector 的值作为新列名
    values_from = baseMean:Label   # 填充 log2FoldChange 的值
  ) %>% mutate(Label = ifelse((Label_HONG31=="Up") & (Label_AAV6=="Up"), "Both up", 
                            ifelse((Label_HONG31=="Down") & (Label_AAV6=="Down"), "Both down","None"))) %>%
    mutate(Label = ifelse((Label_HONG31!=Label_AAV6) & Label=="None", "Only one change", Label))


wide_data%>%as.data.frame%>%group_by( Label_AAV6, Label_HONG31,Label)%>%summarise(count=n())
#add label
label_deg_resdata_filtered_AAV6 <- wide_data %>% 
  #filter(!str_detect(Row.names,"^TR[ABCD]V")) %>% 
  filter(Label != "None") %>% 
  group_by(Label) %>% slice_max(n=5, order_by = abs(log2FoldChange_AAV6)) %>% ungroup 

label_deg_resdata_filtered_HONG31 <- wide_data %>% 
  #filter(!str_detect(Row.names,"^TR[ABCD]V")) %>% 
  filter(Label != "None" & Label != "Only one change" ) %>% 
  group_by(Label) %>% slice_max(n=5, order_by = abs(log2FoldChange_HONG31)) %>% ungroup 

label_deg_resdata_filtered <- rbind(label_deg_resdata_filtered_AAV6,label_deg_resdata_filtered_HONG31) %>% unique

label_deg_resdata_filtered %>% as.data.frame

#cor
wide_data %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31))

cor_time <- signif(cor(wide_data %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
                         wide_data%>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31)),3);cor_time
cor.test(wide_data %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
         wide_data %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31))$p.value

p1 <- ggplot(wide_data , aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31))+
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.6) +
  # 不显著的基因为黑色
  geom_point(data = wide_data  %>% filter(Label=="None"), shape = 21, color = "black", 
             alpha = 0.1, size = 0.1, stroke = 1) 
my_col <- c("Both up" = "#CB181D", "Both down" ="#2927C4", "Only one change" = "#C6A58B")
p2 <- p1 + 
  geom_point(data = wide_data  %>% filter(Label!="None"), 
             aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31, fill = Label, 
                 colour = Label, size = -log10(padj_HONG31), 
                 alpha = -log10(padj_AAV6)),
             shape = 21, stroke = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, col = "orange", alpha=0.5,size=0.7,linetype="dashed")+
  scale_size_continuous(name="-log10(padj_AAV.APR31)",range = c(1.5, 4.5)) +
  scale_alpha_continuous(range = c(0.2, 0.85)) +
  scale_fill_manual(values = my_col ) +
  scale_color_manual(values = my_col )+labs(title=paste0(paste("The correlation of DEGs, no discordant DEG\nPearson r = ",cor_time),
                                                               ", p-value < 2.2e-16 (DEGs only)"))

p_FC <- p2 + 
  ggrepel::geom_text_repel(data= label_deg_resdata_filtered ,
                  aes(x = log2FoldChange_AAV6, y =log2FoldChange_HONG31,label = Row.names),
                  size = 5, color = "black", fontface = 'italic', max.overlaps = 40) +
  #xlim(c(-3, 3)) + 
  #ylim(c(-5, 5)) + 
  xlab("AAV6 D7 vs D16 (log2FC)") +
  ylab("AAV.APR31 D7 vs D16 (log2FC)") +
  theme_bw(base_size = 15) 

ggsave("/home/chengennong/ylp/aav/rnaseq/plot_time/p_FC.pdf",p_FC,height = 5, width = 7.5)



#go-kegg ####

library(clusterProfiler)
library(org.Hs.eg.db)
#top n genes for enrichment
n_genes=100#100


up_genes <- wide_data %>% filter(Label=="Both up")  # 用实际的上调基因替换
up_genes_id <- bitr(up_genes$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db')
down_genes <- wide_data %>% filter(Label=="Both down")  # 用实际的下调基因替换
down_genes_id <- bitr(down_genes$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db')
change_genes <- wide_data %>% filter(Label=="Only one change")  # 用实际的下调基因替换
change_genes_id <- bitr(change_genes$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db')


# GO 富集分析（上调基因）
ego_bp_up <- enrichGO(gene = up_genes_id$SYMBOL,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# GO 富集分析（下调基因）
ego_bp_down <- enrichGO(gene = down_genes_id$SYMBOL,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
# GO 富集分析（上调基因）
ego_mf_up <- enrichGO(gene = up_genes_id$SYMBOL,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# GO 富集分析（下调基因）
ego_mf_down <- enrichGO(gene = down_genes_id$SYMBOL,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "MF",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
# GO 富集分析（上调基因）
ego_cc_up <- enrichGO(gene = up_genes_id$SYMBOL,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# GO 富集分析（下调基因）
ego_cc_down <- enrichGO(gene = down_genes_id$SYMBOL,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

# KEGG 富集分析（上调基因）
kegg_up <- enrichKEGG(gene = up_genes_id$ENTREZID,
                      organism = "hsa",  # 如果是人类数据
                      keyType = "kegg",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# KEGG 富集分析（下调基因）
kegg_down <- enrichKEGG(gene = down_genes_id$ENTREZID,
                        organism = "hsa",  # 如果是人类数据
                        keyType = "kegg",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

pdf("/home/chengennong/ylp/aav/rnaseq/plot_time/go-kegg.pdf", width = 7, height = 11)
dotplot(ego_bp_up, showCategory = 20) + ggtitle("GO BP Enrichment Analysis of Upregulated Genes")
dotplot(ego_bp_down, showCategory = 20) + ggtitle("GO BP Enrichment Analysis of Downregulated Genes")
dotplot(ego_mf_up, showCategory = 20) + ggtitle("GO MF Enrichment Analysis of Upregulated Genes")
dotplot(ego_mf_down, showCategory = 20) + ggtitle("GO MF Enrichment Analysis of Downregulated Genes")
dotplot(ego_cc_up, showCategory = 20) + ggtitle("GO CC Enrichment Analysis of Upregulated Genes")
dotplot(ego_cc_down, showCategory = 20) + ggtitle("GO CC Enrichment Analysis of Downregulated Genes")
dotplot(kegg_up, showCategory = 20) + ggtitle("KEGG Enrichment Analysis of Upregulated Genes")
dotplot(kegg_down, showCategory = 20) + ggtitle("KEGG Enrichment Analysis of Downregulated Genes")
dev.off()


Enrichment_df <- rbind(ego_bp_up@result %>% mutate(Type = "Up_GOBP"),
    ego_bp_down@result %>% mutate(Type = "Down_GOBP"),
    ego_mf_up@result %>% mutate(Type = "Up_GOMF"),
    ego_mf_down@result %>% mutate(Type = "Down_GOMF"),
    ego_cc_up@result %>% mutate(Type = "Up_GOCC"),
    ego_cc_down@result %>% mutate(Type = "Down_GOCC"),
    kegg_up@result %>% mutate(Type = "Up_KEGG"),
    kegg_down@result %>% mutate(Type = "Down_KEGG"))

write.csv(Enrichment_df,"/home/chengennong/ylp/aav/rnaseq/plot_time/Enrichment.csv")

#heatmap ####
counts_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt", header=T) 
length_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/counts.txt", header=T) %>% dplyr::select(1,6)
source("/home/chengennong/ylp/others/scripts/rnaseq.R")
tpm_value <- counts_to_tpm(counts_df0, length_df0)
colnames(tpm_value) <- c(meta_df$Sample)
dim(tpm_value)

cols <- colorRampPalette(colors = c("#8DD3C7", "white", "#C54B8C"))(50)
mat  <- (log(tpm_value+1))[c(label_deg_resdata_filtered$Row.names,"LDHA","KIR2DL3","PDCD1","HAVCR2","LAG3","TOX","ETV1","NR5A2","NR4A2","S1PR1","TCF7","SELL","IL7R","FOXO1","GZMB","GNLY","PRF1"),]
#mat  <- assay(vst_data)[c("LDHA","PDCD1","HAVCR2","LAG3","TOX","NR5A2","S1PR1","TCF7","SELL","IL7R","FOXO1"),]

mat <- mat[,rownames(meta_df %>% filter(str_detect(Group,"elec Cas9 ")))]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))
anno <- as.data.frame((meta_df %>% filter(str_detect(Group,"elec Cas9 ")))[, c("Cell_collection_timepoint","Donor","Group")])
pheatmap::pheatmap(t(mat), annotation_col = anno,
         col = cols,
         width = 12,
         cluster_rows = F,
         cluster_cols = T,
         height = 7,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_time/p_heatmap_gene.pdf")


#Distr
mu <- plyr::ddply(deg_result %>% filter(padj < 0.05), c("Vector"), summarise, grp.mean=mean(log2FoldChange))
p_dist <-ggplot(deg_result %>% filter(padj < 0.05), aes(x=log2FoldChange, color=Vector, fill=Vector)) +
    geom_histogram(position="identity",alpha=0.5)+
    #geom_density(alpha=0.6)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Vector),
               linetype="dashed")+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    labs(title="Log2FC histogram plot",x="Log2FC", y = "Density")+
    #facet_wrap(vars(Time), labeller = "label_both", scales = "free_y", nrow = 1, strip.position = "top") + 
    theme_bw()
library(ggpubr)
#p_dist <- ggdensity(deg_result %>% filter(padj < 0.05), x = "log2FoldChange",
#   add = "mean", rug = TRUE, facet = "Time", ncol = 2, #ylim = c(0, 1),
#   color = "Vector", palette = c("#00AFBB", "#E7B800"))
ggsave("/home/chengennong/ylp/aav/rnaseq/plot_time/p_dist.pdf",p_dist,height = 2, width = 3.5)



#Jaccard
up_gene_df <- deg_result %>% filter(padj < 0.05 & Label=="Up") %>%
  group_by( Vector) #%>% slice_max(n=100, order_by = (log2FoldChange))

down_gene_df <- deg_result %>% filter(padj < 0.05 & Label=="Down") %>%
  group_by( Vector) #%>% slice_max(n=100, order_by = -(log2FoldChange))

jaccard_gene_df <- rbind(up_gene_df %>% mutate(Label="Up"), down_gene_df %>% mutate(Label="Down"))
jaccard_gene_df %>% 
    group_by(Label,  Vector) %>% summarise(count=n()) %>% as.data.frame
colnames(jaccard_gene_df)[1] <- "Gene"
calculate_jaccard <- function(data) {
  # 按 Vector 划分 Gene
  genes_v1 <- data %>% filter(Vector == "AAV6") %>% pull(Gene)
  genes_v2 <- data %>% filter(Vector == "HONG31") %>% pull(Gene)
  # 计算交集和并集
  intersection <- length(intersect(genes_v1, genes_v2))
  union <- length(union(genes_v1, genes_v2))
  # 返回 Jaccard 指数
  if (union == 0) {
    return(NA)  # 避免分母为 0
  } else {
    return(intersection / union)
  }
}
# 按 Label 分组并计算 Jaccard 指数
jaccard_result <- jaccard_gene_df %>%
  group_by(Label) %>%
  summarise(JaccardIndex = calculate_jaccard(cur_data()), .groups = "drop")
jaccard_result


#Num
deg_num_df <- deg_result %>% 
    filter(Label!="None") %>% 
    group_by(Label,  Vector) %>% summarise(count=n()) %>% as.data.frame
deg_num_df <- inner_join(deg_num_df,jaccard_result) %>% mutate(Group=paste(Label));deg_num_df



#p_num <- ggplot(deg_num_df, aes(x= y=count,linetype=Label, fill=paste(Vector))) +
#  geom_bar(stat="identity", color="black", position=position_dodge())+
#  geom_text(aes(y=count, label=count),position = position_dodge(0.9), vjust=1.6, 
#            color="black", size=3.5)+
#  geom_text(aes(y=count+20, label=round(JaccardIndex,3)),position = position_dodge(0.9), vjust=1.6, 
#            color="black", size=3.5)+
#  scale_fill_manual(values=c('#999999','#E69F00'))+
#  #geom_line(aes(linetype=Label))+
#  #facet_wrap(vars(Label), labeller = "label_both", scales = "free_y", nrow = 2, strip.position = "top") + 
#  #facet_wrap(vars(Label, Donor), labeller = "label_both", scales = "free_y", nrow = 2, strip.position = "top") + 
#  theme(legend.position="top") + 
#    theme_bw()

p_num <-ggplot(deg_num_df, aes(x = count, y = Group)) +
  geom_line() +
  geom_point(aes(color = Vector), size = 3) +
  geom_text(aes(x=count+1, label=round(JaccardIndex,3)),position = position_dodge(0.9), vjust=1.6, 
            color="black", size=3.5)+
  scale_color_manual(values=c('#999999','#E69F00'))+ 
  theme_bw()+
  theme(legend.position="top") 
ggsave("/home/chengennong/ylp/aav/rnaseq/plot_time/p_num.pdf",p_num,height = 3, width = 4)
