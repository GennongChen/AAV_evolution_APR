#/home/chengennong/tools/mambaforge/envs/genome/bin/R
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

test_dir <- "/home/chengennong/ylp/aav/screen/screen_out/batch4/mageck_count"

files <- c("Broadgpp_P4_P7.sgrna_summary.txt")

markers <- c("GFP")
gene_sh <- c("CD7","TM9SF2","VPS52","KIAA0319L","NR4A1","KMT2D","HELQ","KPNA6","B2M","GPR108","RRM1","NDNL2","NTC")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"#ColPaired[8]
p_1b_list = list()
for (n_file in 1:length(files)){
    marker <- markers[n_file]
    print(marker)
    a_marker_raw_tests_file <- paste0(test_dir,"/",files[n_file])
    a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
    a_marker_sgrna_tests_table <- a_marker_raw_tests_table %>% 
        mutate(marker=marker) %>% 
        #filter(str_detect(Gene, "^KIAA0319L$|^CD7$|^TM9SF2$|^NR4A1$|^B2M$|^GPR108$|^NTC$")) %>% 
        separate("sgrna",c("sgrna_gene","sgrna_seq"),sep="_") %>% 
        mutate(
        type = case_when(
          #str_detect("IL2|IFNG",Gene) ~ "Control",
          str_detect("^KIAA0319L$|^CD7$|^TM9SF2$|^VPS52$|^NR4A1$|^KMT2D$|^HELQ$|^GPR108$|^KPNA6$",Gene)        ~ "Positive",
          str_detect("^RRM1$|^NDNL2$",Gene)        ~ "Negative",
          str_detect("^NTC$",Gene)        ~ "NTC",
          TRUE                      ~ "other"
        ) 
      )   %>% mutate(sgrna_gene=factor(sgrna_gene,levels=gene_sh)) %>% 
        mutate(type=factor(type,levels=c("NTC","Control","Positive","Negative")))
    a_marker_sgrna_lfc_mean_table <- a_marker_sgrna_tests_table %>% 
        group_by(sgrna_gene,sgrna_seq,type,marker) %>% 
        summarise(LFC=mean(LFC))
    flowcy <- paste0("log2FoldChange(",marker,"_P4","/",marker,"_P7)")
    p_1b_density <-ggplot(a_marker_sgrna_lfc_mean_table, aes(LFC)) +
        geom_density(adjust=1.5, alpha=0.9,fill="grey") +
        scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
        theme_bw()+
        theme(legend.position='none',
              panel.grid =element_blank())+
        labs(x = flowcy,y = "") 
    flowcy <- paste0("log2FoldChange(",marker,"_P4","/",marker,"_P7)")
    p_1b_vline <-ggplot(a_marker_sgrna_lfc_mean_table, aes(LFC, colour = type)) +
        geom_vline(aes(xintercept = LFC, colour = type),size=1, alpha=0.6) +
        facet_grid(sgrna_gene ~ .)+
        scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
        scale_color_manual(values = c("grey","red","blue")) +
        theme_bw()+
        theme(legend.position='none',
              panel.grid =element_blank(),
              strip.background = element_rect(color="white", fill="white", size=1.5))+
        labs(x = flowcy,y = "") 
    
    p_1b <- ggarrange(p_1b_density,p_1b_vline,nrow=2, heights = c(1,2.3))
    p_1b_list[[marker]] <- p_1b
}


#p_1b_total <- ggarrange(plotlist=p_1b_list,ncol=2)
p_1b_total <- ggarrange(plotlist=p_1b_list,ncol=1)
ggsave(p_1b_total,filename="/home/chengennong/ylp/aav/screen/screen_out/batch4/plot/1b.pdf",width=6,height=10)

#print("/home/chengennong/ylp/aav/screen/screen_out/batch4/plot/1b.pdf")
