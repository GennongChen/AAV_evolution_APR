#!/home/chengennong/tools/mambaforge/envs/seurat4/bin/R

#packageVersion("factoextra")
library(tidyverse)

#meta data ####
#Donors=c("Donor 300F", "Donor 1100F")
Time_ps=c("7 Days after elec", "16hours after elec")
Vectors=c("HONG31", "AAV6")

#Donor_i=Donors[1]
#Time_p_i=Time_ps[1]
#Vector_i=Vectors[1]

test_deg_in <- function(Time_p_i,Vector_i) {
    #load data ####
    counts_df <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt", header=T,check.names = F) %>% column_to_rownames("Geneid")
    colnames(counts_df)
    meta_df <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% dplyr::select(Sample,Group,Donor,Cell_collection_timepoint)
    counts_df <- counts_df #%>% dplyr::select(matches(paste(meta_df$Sample, collapse = "|")))
    colnames_count <- sub("/.*", "", colnames(counts_df))
    if (!setequal(colnames_count, meta_df$Sample)) {
      stop("not equal to meta_df")
    }
    colnames(counts_df) <- colnames_count
    counts_df <- counts_df[,meta_df$Sample]
    rownames(meta_df) <- meta_df$Sample
    #names(counts_df)
    head(counts_df)
    #deseq2
    library(DESeq2)  #加载包
    countData <- counts_df %>% as.matrix
    #colData <- data.frame(row.names=meta_df$Sample, Condition=meta_df$Condition)
    colData <- meta_df #%>% select(Group,Cell_collection_timepoint)
    meta_df2 <- read.table("/home/chengennong/ylp/aav/rnaseq/metatable.txt", header=T, sep = "\t") %>% dplyr::select(Sample,Group,Donor,Cell_collection_timepoint) %>% 
      filter(str_detect(Cell_collection_timepoint,Time_p_i) & (str_detect(Group,Vector_i)|str_detect(Group,"elec Cas9$")))
    #正式构建dds矩阵
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~Donor + Group)
    dds <- estimateSizeFactors(dds)
    Nor_countsThreshold <- 2
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds, normalized=TRUE) >= Nor_countsThreshold) >= smallestGroupSize
    #keep2 <- c("3F-31-7-1","11F-0-7-3")
    dds <- dds[keep,]
    counts_df2 <- counts_df[,meta_df2$Sample]
    nrow(dds)
    head(dds)  #查看一下构建好的矩阵
    dds <- DESeqDataSetFromMatrix(countData=counts_df2, colData=meta_df2, design= ~Donor + Group)
    colData(dds)$Group <- relevel(colData(dds)$Group, ref = "elec Cas9")
    dds <- DESeq(dds)   #对原始dds进行normalize
    vst_data <- vst(dds, blind = FALSE)
    pca_assay <- as.data.frame(t(assay(vst_data)))
    #saveRDS(pca_assay, "/home/chengennong/ylp/aav/rnaseq/plot_tech/vst_pca_assay_tech.rds")
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
    resdata_filtered_f <- resdata %>%
      mutate(Comb=paste(Time_p_i,Vector_i,sep="_"),Time=Time_p_i,Vector=Vector_i)
        return(resdata_filtered_f)
}

deg_result_raw <-data.frame()
for (p2 in Time_ps) {
  for (p3 in Vectors) {
    print(paste(p2, p3,sep="_"))
    deg_result_raw <- bind_rows(deg_result_raw, test_deg_in(p2, p3))
  }
}

lfc_thre=1;qval_thre=0.05
deg_result <- deg_result_raw %>% #filter(group == AL[1] & avgExpr!=0) %>% 
      mutate(Label = ifelse((log2FoldChange) > lfc_thre & padj < qval_thre, "Up", 
                            ifelse((log2FoldChange) < -lfc_thre & padj < qval_thre, "Down","None"))) %>% 
      mutate(Label=ifelse(is.na(Label), "None", Label)) %>% 
      mutate(Label=factor(Label, levels= c("Up", "Down", "None"))) %>% 
      arrange(-abs(log2FoldChange), .by_group=TRUE) %>% 
      dplyr::select(Row.names, baseMean, log2FoldChange,padj,Label,Comb,Vector,Time) 

deg_result %>% group_by(Label,Comb,Vector,Time) %>% summarise(count=n())

deg_result %>% filter(!is.na(Label)) %>% head
deg_result%>%filter(Row.names=="WT1")%>%as.data.frame


#Distr
mu <- plyr::ddply(deg_result %>% filter(padj < 0.05), c("Vector","Time"), summarise, grp.mean=mean(log2FoldChange))
p_dist <-ggplot(deg_result %>% filter(padj < 0.05), aes(x=log2FoldChange, color=Vector, fill=Vector)) +
    geom_histogram(position="identity",alpha=0.5)+
    #geom_density(alpha=0.6)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Vector),
               linetype="dashed")+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), labels = c("AAV6", "AAV.APR31"))+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), labels = c("AAV6", "AAV.APR31"))+
    labs(title="DEG distribution",x="log2FoldChange", y = "DEG count")+
    facet_wrap(vars(Time), labeller = "label_both", scales = "free_y", nrow = 1, strip.position = "top") + 
    theme_bw()+ scale_color_manual(name="Condition",values=c("#999999", "#E69F00" ), labels = c("AAV6", "AAV.APR31"))
library(ggpubr)
#p_dist <- ggdensity(deg_result %>% filter(padj < 0.05), x = "log2FoldChange",
#   add = "mean", rug = TRUE, facet = "Time", ncol = 2, #ylim = c(0, 1),
#   color = "Vector", palette = c("#00AFBB", "#E7B800"))
ggsave("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_dist.pdf",p_dist,height = 2.5, width = 7)



deg_result %>% filter(padj < 0.05 & Comb=="16hours after elec_AAV6")
#Jaccard
up_gene_df <- deg_result %>% filter(padj < 0.05 & Label=="Up") %>%
  group_by(Time, Vector) #%>% slice_max(n=100, order_by = (log2FoldChange))
up_gene_df %>% group_by(Comb) %>% summarise(count=n())

down_gene_df <- deg_result %>% filter(padj < 0.05 & Label=="Down") %>%
  group_by(Time, Vector) #%>% slice_max(n=100, order_by = -(log2FoldChange))
down_gene_df %>% group_by(Comb) %>% summarise(count=n())

jaccard_gene_df <- rbind(up_gene_df %>% mutate(Label="Up"), down_gene_df %>% mutate(Label="Down"))
jaccard_gene_df %>% 
    group_by(Label, Time, Vector) %>% summarise(count=n()) %>% as.data.frame
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
# 按 Label, Time 分组并计算 Jaccard 指数
jaccard_result <- jaccard_gene_df %>%
  group_by(Label, Time) %>%
  summarise(JaccardIndex = calculate_jaccard(cur_data()), .groups = "drop")
jaccard_result


#Num
deg_num_df <- deg_result %>% 
    filter(Label!="None") %>% 
    group_by(Label, Time, Vector) %>% summarise(count=n()) %>% as.data.frame
deg_num_df <- inner_join(deg_num_df,jaccard_result) %>% mutate(Group=paste(Time,Label));deg_num_df



p_num <- ggplot(deg_num_df, aes(x=Time, y=count,linetype=Label, fill=paste(Vector))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_text(aes(y=count, label=count),position = position_dodge(0.9), vjust=1.6, 
            color="black", size=3.5)+
  geom_text(aes(y=count+20, label=round(JaccardIndex,3)),position = position_dodge(0.9), vjust=1.6, 
            color="black", size=3.5)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  #geom_line(aes(linetype=Label))+
  #facet_wrap(vars(Label), labeller = "label_both", scales = "free_y", nrow = 2, strip.position = "top") + 
  #facet_wrap(vars(Label, Donor), labeller = "label_both", scales = "free_y", nrow = 2, strip.position = "top") + 
  theme(legend.position="top") + 
    theme_bw()

p_num <-ggplot(deg_num_df, aes(x = count, y = Group)) +
  geom_line() +
  geom_point(aes(color = Vector), size = 3) +
  geom_text(aes(x=count+1, label=round(JaccardIndex,3)),position = position_dodge(0.9), vjust=1.6, 
            color="black", size=3.5)+
  scale_color_manual(values=c('#999999','#E69F00'))+
  theme(legend.position="top") + 
  theme_bw()
ggsave("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_num.pdf",p_num,height = 3, width = 5.5)

#Volcano
resdata_filtered_f <- deg_result %>% mutate(padj=if_else(padj < 1e-10, 1e-10,padj)) %>% 
  #filter(abs(logFC) > 0.05) %>%
  #filter(!str_detect(feature,"^MT-|^RPL|^RPS")) %>% 
  arrange(-abs(log2FoldChange), .by_group=TRUE)
label_deg_resdata_filtered <- resdata_filtered_f %>% 
  #filter(!str_detect(Row.names,"^TR[ABCD]V")) %>% 
  filter(Label != "None") %>% 
  group_by(Label) %>% slice_max(n=10, order_by = abs(log2FoldChange))

p_volcano <- ggplot(resdata_filtered_f, aes(log2FoldChange, -1*log10(padj))) + 
  geom_point(aes(color=Label), alpha=0.7, size=0.5) + 
  ggrepel::geom_text_repel(data = label_deg_resdata_filtered, aes(label = Row.names), size = 2,   max.overlaps = 500) + 
  scale_color_manual(values =c("#CD534CFF","#4C4CFF","grey")) + 
  xlim(c(-10, 10)) + ylim(c(0, 10))+ 
  geom_vline(xintercept=c(-1,1), lty=4, col="black", lwd=0.6) +
  geom_hline(yintercept = -1*log10(0.05), lty=4, col="black", lwd=0.6) + 
  labs(x=expression(log[2](Fold_Change)), y=expression(-log[10](p_adjusted_value)))+
  theme_bw()+
  facet_wrap(vars(Vector,Time))+
  theme(legend.position="top")+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))
ggsave("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_volcano.pdf",p_volcano,height = 6, width = 5)


#Cor
deg_result%>%filter(Row.names=="IGHA2")%>%as.data.frame
wide_data <- deg_result %>% dplyr::select(-Comb) %>%
  pivot_wider(
    names_from = Vector,           # 将 Vector 的值作为新列名
    values_from = baseMean:Label   # 填充 log2FoldChange 的值
  ) %>% mutate(Label = ifelse((Label_HONG31=="Up") & (Label_AAV6=="Up"), "Both up", 
                            ifelse((Label_HONG31=="Down") & (Label_AAV6=="Down"), "Both down","None"))) %>%
    mutate(Label = ifelse((Label_HONG31=="Up" & Label_AAV6=="Down"), "Reverse direction", Label)) %>%
    mutate(Label = ifelse((Label_HONG31=="Down" & Label_AAV6=="Up"), "Reverse direction", Label))


deg_result%>%as.data.frame%>%head
wide_data %>% filter(Label_HONG31=="Up" & Label_AAV6=="Down")
wide_data %>% filter(Label_HONG31=="Down" & Label_AAV6=="Up")
#add label
label_deg_resdata_filtered_AAV6 <- wide_data %>% 
  #filter(!str_detect(Row.names,"^TR[ABCD]V")) %>% 
  filter(Label != "None") %>% 
  group_by(Label, Time) %>% slice_max(n=5, order_by = abs(log2FoldChange_AAV6)) %>% ungroup 

label_deg_resdata_filtered_HONG31 <- wide_data %>% 
  #filter(!str_detect(Row.names,"^TR[ABCD]V")) %>% 
  filter(Label != "None" & Label != "Reverse direction" ) %>% 
  group_by(Label, Time) %>% slice_max(n=5, order_by = abs(log2FoldChange_HONG31)) %>% ungroup 

label_deg_resdata_filtered <- rbind(label_deg_resdata_filtered_AAV6,label_deg_resdata_filtered_HONG31) %>% unique

label_deg_resdata_filtered %>% as.data.frame

#Cor H16
wide_data %>% filter(Time=="16hours after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31))

cor_h16 <- signif(cor(wide_data %>% filter(Time=="16hours after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
                         wide_data %>% filter(Time=="16hours after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31)),3);cor_h16
cor.test(wide_data %>% filter(Time=="16hours after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
         wide_data %>% filter(Time=="16hours after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31))$p.value

#H16
p1 <- ggplot(wide_data %>% filter(Time=="16hours after elec"), aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31))+
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.6) +
  # 不显著的基因为黑色
  geom_point(data = wide_data %>% filter(Time=="16hours after elec") %>% filter(Label=="None"), shape = 21, color = "black", 
             alpha = 0.1, size = 0.1, stroke = 1) 
my_col <- c("Both up" = "#CB181D", "Both down" ="#2927C4", "Reverse direction" = "#C6A58B")
p2 <- p1 + 
  geom_point(data = wide_data  %>% filter(Time=="16hours after elec")%>% filter(Label!="None"), 
             aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31, fill = Label, 
                 colour = Label, size = -log10(padj_HONG31), 
                 alpha = -log10(padj_AAV6)),
             shape = 21, stroke = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, col = "orange", alpha=0.5,size=0.7,linetype="dashed")+
  scale_size_continuous(name="-log10(padj_AAV.APR31)",range = c(1.5, 4.5)) +
  scale_alpha_continuous(range = c(0.2, 0.85)) +
  scale_fill_manual(values = my_col ) +
  scale_color_manual(values = my_col )+labs(title=paste0(paste("The correlation of DEGs\nPearson r = ",cor_h16),
                                                               ", p-value < 2.2e-16 (DEGs only)"))
p_FC <- p2 + 
  ggrepel::geom_text_repel(data= label_deg_resdata_filtered %>% filter(Time=="16hours after elec"),
                  aes(x = log2FoldChange_AAV6, y =log2FoldChange_HONG31,label = Row.names),
                  size = 5, color = "black", fontface = 'italic', max.overlaps = 40) +
  #xlim(c(-3, 3)) + 
  #ylim(c(-5, 5)) + 
  xlab("AAV6 vs Cas9 (log2FC)") +
  ylab("AAV.APR31 vs Cas9 (log2FC)") +
  theme_bw(base_size = 15) 

ggsave("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_FC_H16.pdf",p_FC,height = 5, width = 7)


#Cor D7
cor_d7 <- signif(cor(wide_data %>% filter(Time=="7 Days after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
                         wide_data %>% filter(Time=="7 Days after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31)),3);cor_d7
cor.test(wide_data %>% filter(Time=="7 Days after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_AAV6),
         wide_data %>% filter(Time=="7 Days after elec") %>% filter(!is.na(log2FoldChange_AAV6) & !is.na(log2FoldChange_HONG31) & Label!="None") %>% pull(log2FoldChange_HONG31))$p.value

#D7
p1 <- ggplot(wide_data %>% filter(Time=="7 Days after elec"), aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31))+
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.6) +
  # 不显著的基因为黑色
  geom_point(data = wide_data %>% filter(Time=="7 Days after elec") %>% filter(Label=="None"), shape = 21, color = "black", 
             alpha = 0.1, size = 0.2, stroke = 1) 
my_col <- c("Both up" = "#CB181D", "Both down" ="#2927C4", "Reverse direction" = "#C6A58B")
p2 <- p1 + 
  geom_point(data = wide_data  %>% filter(Time=="7 Days after elec")%>% filter(Label!="None"), 
             aes(x=log2FoldChange_AAV6, y=log2FoldChange_HONG31, fill = Label, 
                 colour = Label, size = -log10(padj_HONG31), 
                 alpha = -log10(padj_AAV6)),
             shape = 21, stroke = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, col = "orange", alpha=0.5,size=0.7,linetype="dashed")+
  scale_size_continuous(name="-log10(padj_AAV.APR31)",range = c(1.5, 4.5)) +
  scale_alpha_continuous(range = c(0.2, 0.85)) +
  scale_fill_manual(values = my_col ) +
  scale_color_manual(values = my_col )+labs(title=paste0(paste("The correlation of DEGs\nPearson r = ",cor_d7),
                                                               "\np-value < 2.2e-16 (DEGs only)"))

p_FC <- p2 + 
  ggrepel::geom_text_repel(data= label_deg_resdata_filtered %>% filter(Time=="7 Days after elec"),
                  aes(x = log2FoldChange_AAV6, y =log2FoldChange_HONG31,label = Row.names),
                  size = 5, color = "black", fontface = 'italic', max.overlaps = 40) +
  #xlim(c(-3, 3)) + 
  #ylim(c(-5, 5)) + 
  xlab("AAV6 vs Cas9 (log2FC)") +
  ylab("AAV.APR31 vs Cas9 (log2FC)") +
  theme_bw(base_size = 15) 

ggsave("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_FC_D7.pdf",p_FC,height = 5, width = 7)

#heatmap ####
counts_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/count1.txt", header=T) 
length_df0 <- read.table("/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata/counts.txt", header=T) %>% dplyr::select(1,6)
source("/home/chengennong/ylp/others/scripts/rnaseq.R")
tpm_value <- counts_to_tpm(counts_df0, length_df0)
colnames(tpm_value) <- c(meta_df$Sample)
dim(tpm_value)

cols <- colorRampPalette(colors = c("#8DD3C7", "white", "#C54B8C"))(50)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))(50)

mat  <- (log(tpm_value+1))[c(#as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names),
  "LDHA","KIR2DL3","PDCD1","HAVCR2","LAG3","TOX","ETV1","NR5A2","NR4A2","S1PR1","TCF7","SELL","IL7R","FOXO1","GZMB","GNLY","PRF1"),]
#mat  <- assay(vst_data)[c("LDHA","PDCD1","HAVCR2","LAG3","TOX","NR5A2","S1PR1","TCF7","SELL","IL7R","FOXO1"),]

for (name in names(modules)) {
  print(modules[[name]])
  module_genes <- modules[[name]] %>% unique
  mat  <- (log(tpm_value+1))[module_genes,]
  anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
  rownames(anno) <- meta_df$Sample
  anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
  mat <- mat[, rownames(anno)]
  mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))
  df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
  result <- df %>%
    group_by(group) %>%
    summarise(across(everything(), mean)) %>%
    select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
    as.matrix()
  rownames(result) <- df$group %>% unique

  anno1 <- anno %>% unique
  rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)
  pheatmap::pheatmap(t(result), annotation_col = anno1,
           col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
           height = 7, width = 8,
           filename = paste0("/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_",name,"_gene.pdf"))
}


modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$Naïve,modules$Cytotoxicity,modules$Exhaustion, modules$`Activation:Effector function`) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]

dim(mat)
as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names)
tpm_value["MT-ATP8",]

anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))


df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)

pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 13, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_all_gene.pdf")


modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$Naïve) %>% unique
#tcell_genes <- c(modules$Naïve,modules$Cytotoxicity,modules$Exhaustion, modules$`Activation:Effector function`) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]
anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))
df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)
pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 2.8, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_naive_gene.pdf")



modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$Cytotoxicity) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]

dim(mat)
as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names)
tpm_value["MT-ATP8",]

anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))


df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)

pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 4.4, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_cyt_gene.pdf")


modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$Exhaustion) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]

dim(mat)
as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names)
tpm_value["MT-ATP8",]

anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))


df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)

pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 3.4, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_ex_gene.pdf")



modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$Glycolysis) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]

dim(mat)
as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names)
tpm_value["MT-ATP8",]

anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))


df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)

pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 3.8, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_la_gene.pdf")

modules <- read_tsv("/home/chengennong/ylp/gbm/metainformation/signature_Chu.txt") %>% 
  as.list %>% map(., function(x) x[!is.na(x)])
tcell_genes <- c(modules$`Activation:Effector function`) %>% unique
#modules$`Oxidative phosphorylation`,modules$`Glycolysis`) %>% unique
mat  <- (log(tpm_value+1))[tcell_genes,]

dim(mat)
as.data.frame(label_deg_resdata_filtered) %>% arrange(log2FoldChange_HONG31)%>% pull(Row.names)
tpm_value["MT-ATP8",]

anno <- as.data.frame((meta_df)[, c("Cell_collection_timepoint","Donor","Group")])
rownames(anno) <- meta_df$Sample
anno <- anno %>% arrange(Cell_collection_timepoint, Donor, Group) 
mat <- mat[, rownames(anno)]
mat <- apply(mat, 1, scales::rescale, to=c(-1, 1))


df <- as.data.frame(cbind(mat,anno)) %>% mutate(group=paste0(Cell_collection_timepoint,Donor,Group))
result <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  select(-Cell_collection_timepoint,-Donor,-Group,-group) %>%
  as.matrix()
rownames(result) <- df$group %>% unique

anno1 <- anno %>% unique
rownames(anno1) <- df$group %>% unique

ann_colors <- list(
  Group = c(`elec Cas9` = "black", `elec Cas9 +AAV6` = '#999999', `elec Cas9 +HONG31 ` = '#E69F00'),  
  Donor = c(`Donor 1100F` = "#86A8E7", `Donor 300F` = "#D16BA5"),
  Cell_collection_timepoint = c(`16hours after elec` = "#40AA40", `7 Days after elec` = "#AAAA40")           
)

pheatmap::pheatmap(t(result), annotation_col = anno1,
         col = cols,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = ann_colors,
         show_colnames = F,
         height = 6.8, width = 6,
         filename = "/home/chengennong/ylp/aav/rnaseq/plot_tech/p_heatmap_act_gene.pdf")
