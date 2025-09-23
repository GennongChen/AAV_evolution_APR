#!/home/chengennong/tools/mambaforge/envs/seurat4/bin/R
library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
library(readxl)
library(ggpubr)
packageVersion("dplyr")


#quantile before or after filtering na?
seq_stat <- function(file){
    id <- str_split(file, pattern = "/",n = 100)[[1]][[8]] #specific path for sample name in metadata!
    donor_name <- meta_info %>% filter(ID==id) %>% pull(Donor)
    round <- meta_info %>% filter(ID==id) %>% pull(Round)
    seqdata <- fread(file, nThread=12) %>% 
        mutate(Count=as.numeric(Count), Donor=donor_name, Round=round, ID=id) %>% 
        filter(!str_detect(Sequence,"\\*"))%>% #filter stop codon
        #filter(case_when(Round=="Round3" ~ Count > 20, #When Round=="Round3", reads number must > 20
        #           T ~ !is.na(Count)))%>%  #Otherwise, filtered
        mutate(Percent = Count/sum(Count))%>% 
        arrange(-Percent) %>% 
        mutate(Rank_Percent=rank(Percent),Rank_Percent_per=rank(Percent)/nrow(.)) %>% 
        #mutate(Labels=ifelse(str_detect(Sequence,".APR..."), "APR_motif", 
        #                  ifelse(Sequence == "GSAQNKD", "AAV6", 
        #                      ifelse(Rank_Percent > max(Rank_Percent)-20, "High", "Low")))
        mutate(Labels=factor(ifelse(str_detect(Sequence,".APR..."), "APR_motif", 
                          ifelse(Sequence == "GSAQNKD", "AAV6", "Non-APR_motif")),
               levels = c("AAV6", "APR_motif", "Non-APR_motif"))) %>%
        mutate(Whe_label=ifelse(rank(-Percent) <= 5 | Labels=="AAV6", "Label", "Non-label"))
    return(seqdata)
}

#read metadata and file path
meta_info <- read_excel("/home/chengennong/ylp/aav/meta_info/20231219_CapIV_metadata.xlsx") %>%
    select(ID, Sample, Round, Donor)
aa_stat_files <- list.files(path="/mnt/cgn/data/aav/aav_batch3_2023_1218/Rawdata", pattern="R1.7mers.aa.stat$", full.names = T, recursive = T)
file_parental="/mnt/cgn/data/aav/aav_batch3_2023_1218/Rawdata/LH007/LH007_R1.7mers.aa.stat"
#read files
parental_aa_stat_df <- seq_stat(file_parental)
aa_stat_df <- map_dfr(aa_stat_files, seq_stat)

#inner_join 3 round of 2 donors with the parental library
all_aa_stat_df <- left_join(aa_stat_df, parental_aa_stat_df %>% select(Sequence, Percent), by="Sequence", suffix = c("_D", "_P")) %>%
    filter(!is.na(Percent_P)) %>%
    mutate(Enrichment = Percent_D/Percent_P)

#plot Donor1 seq
#228B22
#006400
count_cutoff=29; cir_size=15
D1_df <- all_aa_stat_df %>% filter(Donor=="Virus (before selection)"|Donor=="Donor1") %>% arrange(Whe_label, ID)
p_donor1 <- ggplot(D1_df %>% filter(Count > count_cutoff),aes(x = Round, y=log2(Percent_D),size=Enrichment, color=Labels, label = Sequence)) +
        geom_jitter(alpha=0.3, position = position_jitter(seed = 31)) +
        scale_color_manual(values =c("orange", "#006400", "grey")) +
        #scale_color_manual(values =c("blue","orange", "red", "grey")) +
        scale_size(range = c(1, cir_size)) +
        theme_bw() +theme(panel.grid =element_blank())+labs(title="Donor1, Count ≥ 30")
p_donor1_text <- p_donor1 + 
  geom_text_repel(data = subset(D1_df, Whe_label=="Label"|Sequence=="AAPRANE"),
                  aes(label = ifelse(Labels == "AAV6"|Labels == "Non-APR_motif"|Labels == "APR_motif", Sequence, "")),
                  position = position_jitter(seed = 31),
                  size=3, box.padding = 1)
#ggsave(plot=p_donor1, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_d1.png", dpi = 300, width=6,height=5)
ggsave(plot=p_donor1, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_d1.png", dpi = 300, width=9,height=5)
ggsave(plot=p_donor1_text, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_d1_text.png", dpi = 300, width=9,height=5)

#plot Donor3 seq
D3_df <- all_aa_stat_df %>% filter(Donor=="Virus (before selection)"|Donor=="Donor3") %>% arrange(Whe_label, ID)
D3_df %>% filter(log2(Percent_D) > -20)
p_donor3 <- ggplot(D3_df %>% filter(Count > count_cutoff),aes(x = Round, y=log2(Percent_D),size=Enrichment, color=Labels, label = Sequence)) +
        geom_jitter(alpha=0.3, position = position_jitter(seed = 31)) +
        scale_color_manual(values =c("orange", "#006400", "grey")) +
        #scale_color_manual(values =c("blue","orange", "red", "grey")) +
        scale_size(range = c(1, cir_size)) +
        theme_bw()+theme(panel.grid =element_blank())+labs(title="Donor3, Count ≥ 30")
p_donor3_text <- p_donor3 + 
  geom_text_repel(data = subset(D3_df, Whe_label=="Label"|Sequence=="AAPRANE"),
                  aes(label = ifelse(Labels == "AAV6"|Labels == "Non-APR_motif"|Labels == "APR_motif", Sequence, "")),
                  position = position_jitter(seed = 31),
                  size=3, box.padding = 1)
#ggsave(plot=p_donor3, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_d3.png", dpi = 300, width=9,height=5)
ggsave(plot=p_donor3, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_d3.png", dpi = 300, width=9,height=5)
ggsave(plot=p_donor3_text, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_d3_text.png", dpi = 300, width=9,height=5)

#parental all top1%/count>30 per line (APR vs. non-APR) ####
draw_dis <- function(D1_df, donor){
top_seq <- D1_df %>% #group_by(Labels) %>% 
  filter(Round=="Parental" & Count >= count_cutoff) %>% 
  #filter(Round=="Parental" & Rank_Percent_per >= 0.99) %>% 
  pull(Sequence)
data1 <- D1_df %>% filter(Sequence %in% top_seq) %>% 
  group_by(Round,Labels) %>% summarise(Percent_D=sum(Percent_D)) %>% 
  mutate(Round_n=case_when(Round == "Parental" ~ 0,
  Round == "Round1" ~ 1,
  Round == "Round2" ~ 2,
  Round == "Round3" ~ 3))
data1
#area 
p_all_cmt30_a <- ggplot(data1, aes(x=Round_n, y=Percent_D, fill=Labels)) + 
  geom_area()+
  scale_fill_manual(values =c("orange", "#006400", "grey")) +
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Parental all Count >= 30",subtitle = donor)
#line
p_all_cmt30_l <- ggplot(data1, aes(x=Round_n, y=Percent_D, fill=Labels,color=Labels)) +
  geom_line(size=2) +
  geom_point(size=3) +
  scale_color_manual(values =c("orange", "#006400", "grey")) +
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Parental all Count >= 30",subtitle = donor)

#parental APR + non-APR top1000 per line ####
#line
both_top_seq <- D1_df  %>% 
  filter(Round=="Parental") %>% group_by(Labels) %>% 
  slice_max(Rank_Percent_per, n=1e3) %>% 
  pull(Sequence)
both_top_seq %>% length()
data2 <- D1_df %>% filter(Sequence %in% both_top_seq) %>% 
  group_by(Round,Labels) %>% summarise(Percent_D=sum(Percent_D)) %>% 
  mutate(Round_n=case_when(Round == "Parental" ~ 0,
                           Round == "Round1" ~ 1,
                           Round == "Round2" ~ 2,
                           Round == "Round3" ~ 3))
#area

p_top1000_both_a <- ggplot(data2, aes(x=Round_n, y=Percent_D, fill=Labels)) + 
  geom_area()+
  scale_fill_manual(values =c("orange", "#006400", "grey")) +
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Parental APR/non-APR top1000",subtitle = donor)
#line
p_top1000_both_l <- ggplot(data2, aes(x=Round_n, y=Percent_D, fill=Labels,color=Labels)) +
  geom_line(size=2) +
  geom_point(size=3) +
  scale_color_manual(values =c("orange", "#006400", "grey")) +
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Parental APR/non-APR top1000",subtitle = donor)

# enrichment box ####
#test
box_Data <- D1_df %>% filter(Sequence %in% both_top_seq) %>% filter(Round=="Round3") %>% filter(Labels!="AAV6")
wilcox_result <- wilcox.test(Enrichment ~ Labels, data = box_Data)
p_value <- wilcox_result$p.value
if (p_value < 2.2e-16) {p_value <- '< 2.2e-16'};p_value
hline=(D1_df %>%  filter(Round=="Round3") %>% filter(Labels=="AAV6") %>% pull(Enrichment))
p_top1000_enrich_b <-  box_Data %>% 
  mutate(hline=hline) %>% 
  ggplot(data=., aes(x=Labels, y=log1p(Enrichment), fill=Labels,
                     color=Labels)) +
  #geom_violin(inner = NULL,color="black") +
  geom_boxplot(color="black",outlier.alpha=0) +
  geom_jitter(size=0.4, alpha=0.5)+
  scale_fill_manual(values =c("#006400", "grey"))+
  scale_color_manual(values =c("#006400", "grey"))+
  geom_errorbar(aes(y=hline, ymax=hline, ymin=hline), colour="orange",linetype="dashed",width=0.5, alpha=1)+
  #stat_pvalue_manual(data.frame(x = c("APR_motif", "Non-APR_motif"),
  #                                     y = c(max(box_Data$Enrichment) + 1, max(box_Data$Enrichment) + 1),
  #                                     label = paste("p", format(p_value, digits = 3))),
  #                          label = "label", y.position = "y", xmin = "APR_motif", xmax = "Non-APR_motif")+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Parental APR/non-APR top1000",subtitle = donor)
  p_donor <- p_all_cmt30_a/p_all_cmt30_l+p_top1000_both_a+p_top1000_both_l+p_top1000_enrich_b
 return(p_donor)
}
p_d1 <- draw_dis(D1_df,"Donor1")
p_d3 <- draw_dis(D3_df,"Donor3")


wilcox_result
p_do <- ggarrange(p_d1,p_d3,ncol=2,common.legend = T)
ggsave(p_do, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_dis.pdf", width=8,height=13)




#full_join seq stat of 2 donors
com_all_aa_stat_df <- full_join(D1_df %>% filter(Round=="Round3") %>% select(Sequence,Percent_D,Enrichment,Labels),
                                D3_df %>% filter(Round=="Round3") %>% select(Sequence,Percent_D,Enrichment,Labels),
                                by="Sequence",suffix = c("_Donor1", "_Donor3")) %>% 
        mutate_if(is.numeric,coalesce,0)%>% 
        mutate_if(is.factor,coalesce,"Low")%>% 
        mutate(Labels=if_else(Labels_Donor1=="Non-APR_motif"|Labels_Donor3=="Non-APR_motif","Non-APR_motif","NA"))%>% 
        mutate(Labels=if_else(Labels_Donor1=="APR_motif"|Labels_Donor3=="APR_motif","APR_motif",Labels))%>% 
        mutate(Labels=factor(if_else(Labels_Donor1=="AAV6"|Labels_Donor3=="AAV6","AAV6",Labels),
            levels = c("AAV6","APR_motif", "Non-APR_motif"))) %>%
        mutate(Enrichment=sqrt(Enrichment_Donor1*Enrichment_Donor3),Percentage=(Percent_D_Donor1+Percent_D_Donor3)/2) %>% 
  mutate(Whe_label=ifelse(rank(-Percentage) <= 5 | Labels=="AAV6", "Label", "Non-label")) %>% 
  arrange(Whe_label,-Percentage)
aa_stat_df %>% filter(str_detect(Sequence,"LAPRQNE|VVNPAEG|GSAQNKD"))
#options(ggrepel.max.overlaps = 20)

# both seq parental ####
both_top_seq_D1 <- D1_df  %>% 
  filter(Round=="Parental") %>% group_by(Labels) %>% 
  slice_max(Rank_Percent_per, n=1e3) %>% 
  pull(Sequence)
both_top_seq_D3 <- D3_df  %>% 
  filter(Round=="Parental") %>% group_by(Labels) %>% 
  slice_max(Rank_Percent_per, n=1e3) %>% 
  pull(Sequence)
both_top_seq <- c(both_top_seq_D1,both_top_seq_D3) %>% unique
################filtered by percent or FC
cor_donors_per <- signif(cor(com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6") %>% pull(Percent_D_Donor1),
                         com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6") %>% pull(Percent_D_Donor3)),3);cor_donors_per
cor_donors <- signif(cor(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor1),
                         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor3)),3);cor_donors
cor.test(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor1),
         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor3))$p.value


#p_com_per <- ggplot(com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6"),
p_com_per <- ggplot(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6"),
                    #aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),size=Enrichment, color=Labels)) +
                    aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3), color=Labels)) +
  geom_point(alpha=0.3,size=0.5) + 
  geom_text_repel(data=subset(com_all_aa_stat_df, Labels=="AAV6"), aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),label = Sequence), size = 3) +
  #scale_color_manual(values =c("blue","orange", "red", "grey")) +
  scale_color_manual(values =c("orange", "#006400", "grey")) +
  scale_size(range = c(1, 3)) +
  annotate("text", x = -20, y = -10, col = "black", size = 4,
           label = paste("Pearson r = ",cor_donors))+
  geom_smooth(method = "lm", se = TRUE, col = "black", alpha=0.2,size=0.7,linetype="dashed")+
        theme_bw() +
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title=paste0("Donors_cor_percentage, parental APR/non-APR top1000\np-value < 2.2e-16, ",
                                                              paste("Pearson r = ",cor_donors)))#;p_com_per
p_com_per_label <- p_com_per + 
  geom_text_repel(data=subset(com_all_aa_stat_df, Whe_label=="Label"), aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),label = Sequence), size = 3)


#ggsave(p_com_per, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_com_per.pdf", width=8,height=6)
ggsave(p_com_per, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_cor.pdf", width=4.5,height=3.3)
ggsave(p_com_per_label, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_cor_label.pdf", width=4.5,height=3.3)
####
cor_donors_en <- signif(cor(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor1) %>% log1p(),
                         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor3) %>% log1p()),3);cor_donors
cor.test(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor1),
         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor3))

#p_com_en <- ggplot(com_all_aa_stat_df %>% filter(Percentage>10e-5|Labels == "AAV6"),
p_com_en <- ggplot(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6"),
                   aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3), color=Labels)) + 
  geom_point(alpha=0.3,size=0.5) +
  geom_text_repel(data=subset(com_all_aa_stat_df, Labels=="AAV6"), aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3),label = Sequence), size = 3) +
  scale_color_manual(values =c("orange", "#006400", "grey"))+
  theme_bw() +
  annotate("text", x = 2.5, y = 6.5, col = "black", size = 4,
           label = paste("Pearson r = ",cor_donors_en
           ))+
  geom_smooth(method = "lm", se = TRUE, col = "black", alpha=0.2,size=0.7,linetype="dashed")+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title=paste0("Donors_cor_enrichment, parental APR/non-APR top1000\np-value < 2.2e-16, ",
                                                              paste("Pearson r = ",cor_donors_en)))
p_com_en_label <- p_com_en +
    geom_text_repel(data=subset(com_all_aa_stat_df, Whe_label=="Label"), aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3),label = Sequence), size = 3) 


#ggsave(p_com_en, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_com_en.png", dpi = 300, width=8,height=6)
ggsave(p_com_en, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_enrich_cor.pdf", width=4.5,height=3.5)
ggsave(p_com_en_label, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_enrich_cor_label.pdf", width=4.5,height=3.5)

hline_a <- com_all_aa_stat_df %>% filter(Labels == "AAV6") %>% pull(Enrichment) %>% log1p
vline_a <- com_all_aa_stat_df %>% filter(Labels == "AAV6") %>% pull(Percentage) %>% log2
p_com_log2per_en <- ggplot(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6"),
            aes(x = log2(Percentage), y=log1p(Enrichment), color=Labels),linetype = "dashed") +
        geom_point(alpha=0.3) +
  geom_text_repel(data=subset(com_all_aa_stat_df, Labels == "AAV6"), aes(x = log2(Percentage), y=log1p(Enrichment),label = Sequence), size = 3) +
        geom_vline(xintercept = vline_a, lty=4, col="black", lwd=0.6)+
        geom_hline(yintercept = hline_a, lty=4, col="black", lwd=0.6)+
  scale_color_manual(values =c("orange", "#006400", "grey"))+theme_bw()+
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title="Donors_percent_enrichment, parental APR/non-APR top1000",
                                                 subtile = "")
p_com_log2per_en_label <- p_com_log2per_en +
  geom_text_repel(data=subset(com_all_aa_stat_df, Whe_label=="Label"), aes(x = log2(Percentage), y=log1p(Enrichment),label = Sequence), size = 3)

#ggsave(plot=p_com_log2per_en, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_com_log2per_en.png", dpi = 300, width=8,height=6)
ggsave(plot=p_com_log2per_en, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_enrich.pdf", width=5,height=3.5)
ggsave(plot=p_com_log2per_en_label, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_enrich_label.pdf", width=5,height=3.5)

# validatuion_table
candidates_per_df <- (com_all_aa_stat_df) %>%
    arrange(-Percentage)%>%
    slice_head(n = 200)
write.csv(candidates_per_df,file="/home/chengennong/ylp/aav/meta_info/candidates_per_CapIV.csv")

candidates_df <- subset(com_all_aa_stat_df,(Labels == "AAV6"|Enrichment>200)) %>%
    arrange(-Percentage)%>%
    slice_head(n = 20)
write.csv(candidates_df,file="/home/chengennong/ylp/aav/meta_info/candidates_CapIV.csv")

candidates_fcgt200_top50_meanrank_df <- subset(com_all_aa_stat_df,(Enrichment>200)) %>%
    mutate(Rank_Percentage=rank(-Percentage),Rank_Enrichment=rank(-Enrichment),Rank_mean=rank(-Percentage)/2+rank(-Enrichment))%>%
    arrange(Rank_mean)
candidates_fcgt200_top50_meanrank_df %>% filter((str_detect(Sequence,".APR...")))
write.csv(candidates_fcgt200_top50_meanrank_df,file="/home/chengennong/ylp/aav/meta_info/candidates_rank_CapIV_enMt200.csv")

candidates_fcgt200_APRonly_meanrank_df <- subset(com_all_aa_stat_df,(Enrichment>200 & Labels == "APR_motif")) %>%
    mutate(Rank_Percentage=rank(-Percentage),Rank_Enrichment=rank(-Enrichment),Rank_mean=rank(-Percentage)/2+rank(-Enrichment))%>%
    arrange(Rank_mean)
candidates_fcgt200_APRonly_meanrank_df %>% filter((str_detect(Sequence,".APR...")))
write.csv(candidates_fcgt200_APRonly_meanrank_df,file="/home/chengennong/ylp/aav/meta_info/candidates_rank_CapIV_enMt200_APRonly.csv")


candidates_fcgt300_top50_meanrank_df <- subset(com_all_aa_stat_df,(Enrichment>300)) %>%
    mutate(Rank_Percentage=rank(-Percentage),Rank_Enrichment=rank(-Enrichment),Rank_mean=rank(-Percentage)/2+rank(-Enrichment))%>%
    arrange(Rank_mean)
candidates_fcgt300_top50_meanrank_df %>% filter((str_detect(Sequence,".APR...")))
write.csv(candidates_fcgt300_top50_meanrank_df,file="/home/chengennong/ylp/aav/meta_info/candidates_rank_CapIV.csv")
##motif
library(ggseqlogo)

#?ggseqlogo
#FC>300 & Rank_mean 50
ggsave(plot=ggseqlogo(candidates_fcgt300_top50_meanrank_df$Sequence), filename="/home/chengennong/ylp/aav/fig_CapIV/motif_gt300_top50.pdf",width=8,height=4)
candidates_fcgt300_top50_meanrank_df%>%group_by(Labels)%>%summarise(count=n())
#plot protein sequence logo
FC_gt300_all_df <- subset(com_all_aa_stat_df, (Enrichment>300))
ggsave(plot=ggseqlogo(FC_gt300_all_df$Sequence), filename="/home/chengennong/ylp/aav/fig_CapIV/ye/motif_both100_all.pdf",width=8,height=4)
FC_gt300_all_df%>%group_by(Labels)%>%summarise(count=n())

top_seq <- D1_df %>% #group_by(Labels) %>% 
  filter(Round=="Parental" & Count >= count_cutoff) %>% 
  #filter(Round=="Parental" & Rank_Percent_per >= 0.99) %>% 
  pull(Sequence)

#top1000 round3
top_n_var <- 200
both_top_seq_logo <- com_all_aa_stat_df  %>% 
  #group_by(Labels) %>% 
  slice_max(Percentage , n=top_n_var) %>% 
  pull(Sequence)
per_2_donor <- com_all_aa_stat_df  %>% 
  slice_max(Percentage , n=top_n_var) %>% 
  summarise(Percent_D_Donor1=sum(Percent_D_Donor1), Percent_D_Donor3=sum(Percent_D_Donor3))
p_both_top100_motif <- ggseqlogo(both_top_seq_logo)+labs(title=paste0("Round3 all top200\nDonor1: ",signif(per_2_donor[1],3),", Donor2: ",signif(per_2_donor[2],3)))
#ggsave(p_both_top100_motif, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/motif_both_top100_all.pdf",width=8,height=4)
top_seq_D1 <- D1_df %>% #group_by(Labels) %>% 
  filter(Round=="Round3") %>% 
  slice_max(Percent_D, n=top_n_var) %>% 
  pull(Sequence)
top_seq_D3 <- D3_df %>% #group_by(Labels) %>% 
  filter(Round=="Round3") %>% 
  slice_max(Percent_D, n=top_n_var) %>% 
  pull(Sequence)
per_donor1 <- D1_df  %>% 
  slice_max(Percent_D , n=top_n_var) %>% 
  summarise(Percent_D=sum(Percent_D))
per_donor3 <- D3_df  %>% 
  slice_max(Percent_D , n=top_n_var) %>% 
  summarise(Percent_D=sum(Percent_D))
p_d1_top100_motif <- ggseqlogo(top_seq_D3)+labs(title=paste0("Round3 all top200\nDonor1: ",signif(per_donor1[1],3)))
p_d3_top100_motif <- ggseqlogo(top_seq_D1)+labs(title=paste0("Round3 all top200\nDonor3: ",signif(per_donor3[1],3)))
p_nmotif <- ggarrange(p_d1_top100_motif,p_d3_top100_motif,p_both_top100_motif,ncol=3,common.legend = T)
ggsave(p_nmotif, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/p_motif.pdf",width=8,height=2)
#p_den<- ggplot(com_all_aa_stat_df %>%filter((Percentage>0.0001)), aes(x=Enrichment, color=Labels)) +
#  geom_density()+
#  #geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
#  #           linetype="dashed")+
#  #labs(title="Weight density curve",x="Weight(kg)", y = "Density")+ 
#  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#  theme_classic()
#ggsave(plot=p_den, filename="/home/chengennong/ylp/aav/fig_CapIV/density.pdf",width=4,height=2)
#stat
seq_meta_stat <- function(seqdata){
    Unique_reads <- seqdata %>% 
                    group_by(Round, Donor) %>%
                    summarise(Unique_reads = n())
    Total_reads <- seqdata %>% 
                   group_by(Round, Donor) %>%
                   summarise(Total_reads = sum(Count))
    stat_df <- inner_join(Total_reads,Unique_reads,by=c("Round", "Donor"))
    return(stat_df)
}
seq_meta_stat_df <- seq_meta_stat(aa_stat_df)
write.csv(seq_meta_stat_df,file="/home/chengennong/ylp/aav/table/reads_stat_CapIV.csv")
bar_df<-read.csv(file="/home/chengennong/ylp/aav/table/stat_bar_CapIV.csv")
#plot meta stat
p_stat_t <- ggplot(bar_df %>% pivot_longer(ends_with("reads"), names_to = "Reads_type", values_to = "Reads_number")%>% filter(Reads_type=="Total_reads"),
                    aes(x=Donor, y=Reads_number, fill=Round)) +
   geom_bar(stat="identity", position=position_dodge())+ 
   #scale_fill_brewer(palette="Greens") +labs(title="Total reads")+theme_bw()
   scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(12, "Greens"))(8)[3:6])+labs(title="Total reads")+theme_bw()

p_stat_u <- ggplot(bar_df %>% pivot_longer(ends_with("reads"), names_to = "Reads_type", values_to = "Reads_number") %>% filter(Reads_type=="Unique_reads"),
                    aes(x=Donor, y=Reads_number, fill=Round)) + 
   geom_bar(stat="identity", position=position_dodge())+ 
   #scale_fill_brewer(palette="Greens",)+labs(title="Unique reads")+theme_bw()
   scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(12, "Greens"))(8)[3:6])+labs(title="Unique reads")+theme_bw()

#ggsave(plot=(p_stat_t+p_stat_u & theme(legend.position = "right")) +plot_layout(guides = "collect"), filename="/home/chengennong/ylp/aav/fig_CapIV/stat.pdf",width=8,height=3)
ggsave(plot=(p_stat_t+p_stat_u & theme(legend.position = "right")) +plot_layout(guides = "collect"), filename="/home/chengennong/ylp/aav/fig_CapIV/ye/stat.pdf",width=8,height=3)




# both seq r3 ####
both_top_seq_D1 <- D1_df  %>% 
  filter(Round=="Round3") %>% #group_by(Labels) %>% 
  slice_max(Rank_Percent_per, n=1e3) %>% 
  pull(Sequence)
both_top_seq_D3 <- D3_df  %>% 
  filter(Round=="Round3") %>% #group_by(Labels) %>% 
  slice_max(Rank_Percent_per, n=1e3) %>% 
  pull(Sequence)
both_top_seq <- c(both_top_seq_D1,both_top_seq_D3) %>% unique
################filtered by percent or FC
cor_donors_per <- signif(cor(com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6") %>% pull(Percent_D_Donor1),
                         com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6") %>% pull(Percent_D_Donor3)),3);cor_donors_per
cor_donors <- signif(cor(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor1),
                         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor3)),3);cor_donors
cor.test(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor1),
         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Percent_D_Donor3))$p.value


#p_com_per <- ggplot(com_all_aa_stat_df %>% filter((Enrichment)>50|Labels == "AAV6"),
p_com_per <- ggplot(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6"),
                    #aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),size=Enrichment, color=Labels)) +
                    aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3), color=Labels)) +
  geom_point(alpha=0.3,size=0.5) + 
  geom_text_repel(data=subset(com_all_aa_stat_df, Labels=="AAV6"), aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),label = Sequence), size = 3) +
  #scale_color_manual(values =c("blue","orange", "red", "grey")) +
  scale_color_manual(values =c("orange", "#006400", "grey")) +
  scale_size(range = c(1, 3)) +
  annotate("text", x = -20, y = -10, col = "black", size = 4,
           label = paste("Pearson r = ",cor_donors))+
  geom_smooth(method = "lm", se = TRUE, col = "black", alpha=0.2,size=0.7,linetype="dashed")+
        theme_bw() +
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title=paste0("Donors_cor_percentage, parental APR/non-APR top1000\np-value < 2.2e-16, ",
                                                              paste("Pearson r = ",cor_donors)))#;p_com_per
p_com_per_label <- p_com_per + 
  geom_text_repel(data=subset(com_all_aa_stat_df, Whe_label=="Label"), aes(x = log2(Percent_D_Donor1), y=log2(Percent_D_Donor3),label = Sequence), size = 3)


#ggsave(p_com_per, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_com_per.pdf", width=8,height=6)
ggsave(p_com_per, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_cor_r3.pdf", width=4.5,height=3.3)
ggsave(p_com_per_label, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_per_cor_label_r3.pdf", width=4.5,height=3.3)
####
cor_donors_en <- signif(cor(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor1) %>% log1p(),
                         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor3) %>% log1p()),3);cor_donors
cor.test(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor1),
         com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6") %>% pull(Enrichment_Donor3))$p.value

#p_com_en <- ggplot(com_all_aa_stat_df %>% filter(Percentage>10e-5|Labels == "AAV6"),
p_com_en <- ggplot(com_all_aa_stat_df %>% filter(Sequence %in% both_top_seq|Labels == "AAV6"),
                   aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3), color=Labels)) + 
  geom_point(alpha=0.3,size=0.5) +
  geom_text_repel(data=subset(com_all_aa_stat_df, Labels=="AAV6"), aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3),label = Sequence), size = 3) +
  scale_color_manual(values =c("orange", "#006400", "grey"))+
  theme_bw() +
  annotate("text", x = 2.5, y = 6.5, col = "black", size = 4,
           label = paste("Pearson r = ",cor_donors_en
           ))+
  geom_smooth(method = "lm", se = TRUE, col = "black", alpha=0.2,size=0.7,linetype="dashed")+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid =element_blank(),
        axis.line = element_line(size=0.3))+labs(title=paste0("Donors_cor_enrichment, parental APR/non-APR top1000\np-value < 2.2e-16, ",
                                                              paste("Pearson r = ",cor_donors_en)))
p_com_en_label <- p_com_en +
    geom_text_repel(data=subset(com_all_aa_stat_df, Whe_label=="Label"), aes(x = log1p(Enrichment_Donor1), y=log1p(Enrichment_Donor3),label = Sequence), size = 3) 


#ggsave(p_com_en, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_com_en.png", dpi = 300, width=8,height=6)
ggsave(p_com_en, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_enrich_cor_r3.pdf", width=4.5,height=3.5)
ggsave(p_com_en_label, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_enrich_cor_label_r3.pdf", width=4.5,height=3.5)
