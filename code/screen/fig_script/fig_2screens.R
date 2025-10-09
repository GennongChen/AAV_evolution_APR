#!/home/chengennong/tools/mambaforge/envs/seurat4/bin/R
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(colorspace)

test_dir <- c("/home/chengennong/ylp/aav/screen/screen_out/batch4/mageck_count","/home/chengennong/ylp/aav/screen/screen_out/batch5/mageck_count")

files <- c("Broadgpp_P4_P7.sgrna_summary.txt", "Broadgpp_P4_P7.sgrna_summary.txt")
gene_files <- c("Broadgpp_P4_P7.gene_summary.txt","Broadgpp_P4_P7.gene_summary.txt")
markers <- c("Jukrat_GFP_P7","Nalm6_GFP_P7")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"
###Jukrat_GFP_P7 1d
df_list <- list()
marker <- markers[1]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir[1],"/",gene_files[1])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id)
a_marker_raw_tests_file <- paste0(test_dir[1],"/",files[1])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)


a_marker_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
  mutate(marker=marker) %>% 
  separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
  inner_join(a_marker_raw_gene_tests_table) %>% 
  group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc)  %>% 
  mutate(med_LFC=median(LFC)) %>% 
  group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc,neg.lfc,med_LFC)  %>% 
  summarise(mean_donor_LFC=mean(LFC)) %>% 
  mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive Hit",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative Hit",
      TRUE                      ~ "Not a Hit"
    ))   %>% 
  pivot_wider(names_from = sgrna_donor, values_from = mean_donor_LFC) %>% 
  mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
  mutate(type=factor(type,levels=c("Not a Hit","Positive Hit","Negative Hit"))) %>% 
  mutate(mean_donor_LFC=(r0+r1+r2)/3)
#mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))
head(a_marker_sgrna_gene_tests_table) %>% as.data.frame
df_list[[marker]] <- a_marker_sgrna_gene_tests_table 


###Nalm6_GFP_P7 1d
marker <- markers[2]
print(marker)
b_marker_raw_gene_tests_file <- paste0(test_dir[2],"/",gene_files[2])
b_marker_raw_gene_tests_table <- read.table(b_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id)
b_marker_raw_tests_file <- paste0(test_dir[2],"/",files[2])
b_marker_raw_tests_table <- read.table(b_marker_raw_tests_file,sep="\t",header=1)

b_marker_sgrna_gene_tests_table <- b_marker_raw_tests_table %>% 
  mutate(marker=marker) %>% 
  separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
  inner_join(b_marker_raw_gene_tests_table) %>% 
  group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc)  %>% 
  mutate(med_LFC=median(LFC)) %>% 
  group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc,neg.lfc,med_LFC)  %>% 
  summarise(mean_donor_LFC=mean(LFC)) %>% 
  mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive Hit",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative Hit",
      TRUE                      ~ "Not a Hit"
    ))   %>% 
  pivot_wider(names_from = sgrna_donor, values_from = mean_donor_LFC) %>% 
  mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
  mutate(type=factor(type,levels=c("Not a Hit","Positive Hit","Negative Hit"))) %>% 
  mutate(mean_donor_LFC=(r0+r1+r2)/3)
#mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))

df_list[[marker]] <- b_marker_sgrna_gene_tests_table 

# 2 markers
a_2marker_sgrna_gene_tests_fea_table0 <- Reduce(full_join,df_list)  %>% ungroup %>% 
  select(sgrna_gene,type,marker) %>% mutate(type=as.character(type)) %>% 
  pivot_wider(names_from = marker, values_from = type) %>%  
  mutate(type = case_when(
    (str_detect("Not a Hit",Jukrat_GFP_P7) & str_detect("Not a Hit",Nalm6_GFP_P7))  ~ "Not a Hit",
    (str_detect("Positive Hit",Jukrat_GFP_P7) & !str_detect("Positive Hit",Nalm6_GFP_P7))  ~ "Jukrat_GFP_P7 Only",
    (str_detect("Positive Hit",Nalm6_GFP_P7) & !str_detect("Positive Hit",Jukrat_GFP_P7))  ~ "Nalm6_GFP_P7 Only",
    (str_detect("Positive Hit",Jukrat_GFP_P7) & str_detect("Positive Hit",Nalm6_GFP_P7))  ~ "Both Screens",
    TRUE  ~ "Not a Hit")) %>% 
  #mutate(type = case_when(
  #      (str_detect("Not a Hit",Jukrat_GFP_P7) & str_detect("Not a Hit",Nalm6_GFP_P7))  ~ "Not a Hit",
  #      (str_detect("Positive Hit|Negative Hit",Jukrat_GFP_P7) & str_detect("Not a Hit",Nalm6_GFP_P7))  ~ "Jukrat_GFP_P7 Only",
  #      (str_detect("Positive Hit|Negative Hit",Nalm6_GFP_P7) & str_detect("Not a Hit",Jukrat_GFP_P7))  ~ "Nalm6_GFP_P7 Only",
  #      (str_detect("Positive Hit|Negative Hit",Jukrat_GFP_P7) & str_detect("Positive Hit|Negative Hit",Nalm6_GFP_P7))  ~ "Both Screens",
  #      TRUE  ~ "other")) %>% 
  select(sgrna_gene, type)

spe_label <- c("TM9SF2", "GPR108", "VPS52", "VPS51", "VPS41", "VPS16", "ATP2C1", "TRAPPC1", "COG4", "OSBPL8", "SLC35B2")
a_2marker_sgrna_gene_tests_fea_table <- Reduce(full_join,df_list)  %>% ungroup %>% 
  select(sgrna_gene,marker,mean_donor_LFC) %>% 
  inner_join(a_2marker_sgrna_gene_tests_fea_table0) %>% 
  pivot_wider(names_from = marker, values_from = mean_donor_LFC) %>% 
  mutate(type=factor(type,levels=c("Both Screens","Jukrat_GFP_P7 Only","Nalm6_GFP_P7 Only","Not a Hit"))) %>% 
  mutate(label = ifelse(abs(Jukrat_GFP_P7) > 1.8 | abs(Nalm6_GFP_P7) > 2 |sgrna_gene %in% spe_label , sgrna_gene, ""))

a_2marker_sgrna_gene_tests_fea_table0 %>% filter(sgrna_gene=="GPR108")
a_2marker_sgrna_gene_tests_fea_table %>% filter(sgrna_gene=="GPR108")
a_2marker_sgrna_gene_tests_fea_table$type %>% unique

p_1e_total <-ggplot(a_2marker_sgrna_gene_tests_fea_table, aes(Jukrat_GFP_P7, Nalm6_GFP_P7)) +
  geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
  geom_density_2d(color="black") +
  geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
  geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
  geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,max.overlaps = 9000,seed = 7654) +
  #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
  theme_bw()+
  labs(x = "Jukrat-GFP-LFC(P4/P7)",y = "Nalm6-GFP-LFC(P4/P7)", title="CRISPR KO Screens",fill="Screen Hit") +
  scale_fill_manual(values = c("purple","red","blue","grey")) +
  scale_color_manual(values = darken(c("purple","red","blue","grey"), 0.3))+
  theme(legend.position='right',
        panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框

ggsave(p_1e_total,filename="/home/chengennong/ylp/aav/screen/screen_out/fig_script/plot/1e2.pdf",width=7,height=5)


#go-kegg ####

#a_2marker_sgrna_gene_tests_fea_table0 %>% filter(type=="Jukrat_GFP_P7 Only") %>% pull(sgrna_gene) %>% unique
#
#degs_df <-  a_marker_sgrna_gene_tests_table %>%
#  mutate(Label = ifelse((pos.lfc) > 0.5 & pos.fdr < 0.05, "Up", 
#                        ifelse((neg.lfc) < -0.5 & neg.fdr < 0.05, "Down","None")))%>% 
#  group_by(Label) %>% filter(Label!="None") %>% arrange(neg.rank)
#
#library(clusterProfiler)
#library(org.Hs.eg.db)
#
#degs_id_df <- bitr(degs_df$sgrna_gene,'SYMBOL','ENTREZID','org.Hs.eg.db')%>% dplyr::rename(Row.names=SYMBOL) %>% ## SYMBOL2ENTREZID
#  left_join(degs_df%>% dplyr::rename(Row.names=sgrna_gene)) #%>% head
#
#
#count(degs_id_df,Label)
#
#
#up_genes <- degs_id_df %>% filter(Label=="Up")  # 用实际的上调基因替换
#up_genes_id <- bitr(up_genes$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db')
#down_genes <- degs_id_df %>% filter(Label=="Down")  # 用实际的下调基因替换
#down_genes_id <- bitr(down_genes$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db')
#
#library(clusterProfiler)
#library(org.Hs.eg.db)
## GO 富集分析（上调基因）
#ego_bp_up <- enrichGO(gene = up_genes_id$SYMBOL,
#                   OrgDb = org.Hs.eg.db,
#                   keyType = "SYMBOL",
#                   ont = "all",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
#                   pAdjustMethod = "BH",
#                   pvalueCutoff = 0.05,
#                   qvalueCutoff = 0.05)
#
## GO 富集分析（下调基因）
#ego_bp_down <- enrichGO(gene = down_genes_id$SYMBOL,
#                     OrgDb = org.Hs.eg.db,
#                     keyType = "SYMBOL",
#                     ont = "all",  # 选择 Biological Process, 还可以选择 "MF" 或 "CC"
#                     pAdjustMethod = "BH",
#                     pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.05)
#
#
## KEGG 富集分析（上调基因）
#kegg_up <- enrichKEGG(gene = up_genes_id$ENTREZID,
#                      organism = "hsa",  # 如果是人类数据
#                      keyType = "kegg",
#                      pAdjustMethod = "BH",
#                      pvalueCutoff = 0.05,
#                      qvalueCutoff = 0.05)
#
## KEGG 富集分析（下调基因）
#kegg_down <- enrichKEGG(gene = down_genes_id$ENTREZID,
#                        organism = "hsa",  # 如果是人类数据
#                        keyType = "kegg",
#                        pAdjustMethod = "BH",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.05)
#
#pdf("/home/chengennong/ylp/aav/screen/screen_out/batch4/plot/1d.pdf", width = 7, height = 8)
#dotplot(ego_bp_up, showCategory = 20) + ggtitle("GO Enrichment Analysis of Upregulated Genes")
#dotplot(ego_bp_down, showCategory = 20) + ggtitle("GO Enrichment Analysis of Downregulated Genes")
#dotplot(kegg_up, showCategory = 20) + ggtitle("KEGG Enrichment Analysis of Upregulated Genes")
#dotplot(kegg_down, showCategory = 20) + ggtitle("KEGG Enrichment Analysis of Downregulated Genes")
#dev.off()
#
#ego_bp_up@result
#
#Enrichment_df <- rbind(ego_bp_up@result %>% mutate(Type = "Up_GOBP"),
#    ego_bp_down@result %>% mutate(Type = "Down_GOBP"),
#    kegg_up@result %>% mutate(Type = "Up_KEGG"),
#    kegg_down@result %>% mutate(Type = "Down_KEGG"))
#
#write.csv(Enrichment_df,"/home/chengennong/ylp/others/pym/rnaseq/plot/Enrichment.csv")
#