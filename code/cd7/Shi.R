#.libPaths("/home/chengennong/tools/mambaforge/envs/SCP_env/lib/R/library")


library(tidyverse)
library(Seurat)

#pseudo bulk shi####
exp_cd7_aavr <- read.csv("/home/chengennong/ylp/aav/cd7_exp/fig/Shi_exp_df.csv")
##tissue

exp_cd7_aavr_long <- exp_cd7_aavr %>% 
  pivot_longer(cols = CD7:KIAA0319L, names_to = "Gene name", values_to = "nExpression")
exp_cd7_aavr_mean <- exp_cd7_aavr_long %>% group_by(tissue,system,`Gene name`) %>% summarise(nExpression=mean(nExpression))

expr_df <- exp_cd7_aavr_mean %>% filter(`Gene name`=="CD7") %>% arrange(nExpression)
expr_df$tissue %>% unique
expr_df_cd7_aavr <- exp_cd7_aavr_mean %>% filter(`Gene name`=="CD7" | `Gene name`=="KIAA0319L") %>% 
  mutate(Tissue=factor(tissue,levels = (expr_df$tissue %>% unique)))

p_shi_tissue_bar<-ggplot(data=expr_df_cd7_aavr, aes(x=Tissue, y=nExpression)) +
  geom_bar(stat="identity", fill='#C70039')+
  facet_wrap(~`Gene name`)+ coord_flip() +
  theme_minimal()+
  labs(title = "Human pseudo-bulk level (Shi et.al)",xlab="Log2(Nornmalized Count + 1)")
ggsave(filename="/home/chengennong/ylp/aav/cd7_exp/fig/bar_CD7_AAVR_shi.pdf", p_shi_tissue_bar, width =3, height = 6)

p_prop <- ggplot(exp_cd7_aavr %>% mutate(majorCluster=factor(majorCluster,levels = c("ILC","CD8T","CD4T","B","Myeloid","Endothelial","Epithelial","Stromal"))) %>% 
                   mutate(Tissue=factor(tissue,levels = (expr_df$tissue %>% unique))), 
                 aes(x=Tissue, fill=majorCluster)) + geom_bar(position = "fill") +
  scale_fill_manual(values = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f"))+
  RotatedAxis()+ggplot2::coord_flip() +theme_classic()
ggsave(filename="/home/chengennong/ylp/aav/cd7_exp/fig/bar_prop.pdf", p_prop, width =4, height = 5)

##system
exp_cd7_aavr_mean <- exp_cd7_aavr_long %>% group_by(system,`Gene name`) %>% summarise(nExpression=mean(nExpression))
expr_df <- exp_cd7_aavr_mean %>% filter(`Gene name`=="CD7") %>% arrange(nExpression)
expr_df$system %>% unique
expr_df_cd7_aavr <- exp_cd7_aavr_mean %>% filter(`Gene name`=="CD7" | `Gene name`=="KIAA0319L") %>% 
  mutate(System=factor(system,levels = (expr_df$system %>% unique)))

p_shi_sys_bar<-ggplot(data=expr_df_cd7_aavr, aes(x=System, y=nExpression)) +
  geom_bar(stat="identity", fill='#C70039')+
  facet_wrap(~`Gene name`)+ coord_flip() +
  theme_minimal()+
  labs(title = "Human pseudo-bulk system level (Shi et.al)")
ggsave(filename="/home/chengennong/ylp/aav/cd7_exp/fig/bar_CD7_AAVR_shi_sys.pdf", p_shi_sys_bar, width =3, height = 3)

p_prop_sys <- ggplot(exp_cd7_aavr %>% mutate(majorCluster=factor(majorCluster,levels = c("ILC","CD8T","CD4T","B","Myeloid","Endothelial","Epithelial","Stromal"))) %>% 
                       mutate(System=factor(system,levels = (expr_df$system %>% unique))), 
                 aes(x=System, fill=majorCluster)) + geom_bar(position = "fill") +
  scale_fill_manual(values = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f"))+
  RotatedAxis()+ggplot2::coord_flip() +theme_classic()
ggsave(filename="/home/chengennong/ylp/aav/cd7_exp/fig/bar_prop_sys.pdf", p_prop_sys, width =4, height = 2.5)


#HPA ####
hpa_exp_df <- read_tsv("/home/chengennong/ylp/aav/cd7_exp/rna_tissue_hpa.tsv")
hpa_exp_cd7 <- hpa_exp_df %>% filter(`Gene name`=="CD7") %>% arrange(nTPM)
hpa_exp_cd7$Tissue
hpa_exp_cd7_aavr <- hpa_exp_df %>% filter(`Gene name`=="CD7" | `Gene name`=="KIAA0319L") %>% 
  mutate(Tissue=factor(Tissue,levels = hpa_exp_cd7$Tissue))

p_hpa_bar<-ggplot(data=hpa_exp_cd7_aavr, aes(x=Tissue, y=log2(nTPM+1))) +
  geom_bar(stat="identity", fill='#C70039')+
  facet_wrap(~`Gene name`)+ coord_flip() +
  theme_minimal()+
  labs(title = "Human bulk level (HPA)")
ggsave(filename="/home/chengennong/ylp/aav/cd7_exp/fig/bar_CD7_AAVR_hpa.pdf", p_hpa_bar, width =3, height = 5.6)







