#!/usr/bin/R
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(colorspace)

test_dir <- "/home/chengennong/ylp/aav/screen/screen_out/batch4/mageck_count"

files <- c("Broadgpp_P4_P7.sgrna_summary.txt", "Broadgpp_P4_baseline.sgrna_summary.txt")
gene_files <- c("Broadgpp_P4_P7.gene_summary.txt","Broadgpp_P4_baseline.gene_summary.txt")
markers <- c("GFP_P7","GFP_baseline")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"

marker <- markers[1]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[1])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id)
a_marker_raw_tests_file <- paste0(test_dir,"/",files[1])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)


a_marker_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc,pos.p.value,pos.score)  %>% 
    mutate(med_LFC_sg=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,pos.lfc,neg.lfc,med_LFC_sg,pos.p.value,pos.score)  %>% 
    summarise(mean_donor_LFC=mean(LFC)) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      pos.lfc > 0.5 & pos.fdr < 0.05        ~ "Positive Hit",
      neg.lfc < -0.5 & neg.fdr < 0.05        ~ "Negative Hit",
      TRUE                      ~ "Not a Hit"
    ))   %>% 
    pivot_wider(names_from = sgrna_donor, values_from = mean_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Positive Hit","Negative Hit","Not a Hit"))) %>% 
    mutate(mean_donor_LFC=(r0+r1+r2)/2) 
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))
head(a_marker_sgrna_gene_tests_table %>% filter(sgrna_gene=="CD7")) %>% as.data.frame
a_marker_sgrna_gene_tests_table%>% arrange(-pos.lfc)
count_cutoff=0; cir_size=10;options(ggrepel.max.overlaps = Inf)

p_bubble <- ggplot(a_marker_sgrna_gene_tests_table,aes(x = marker, y=-log2(pos.fdr),size=2^(pos.lfc), color=type, label = sgrna_gene)) +
        geom_jitter(alpha=0.3, position = position_jitter(seed = 31)) +
  geom_text_repel(#data = subset(a_marker_sgrna_gene_tests_table, type=="Positive Hit"),
                  aes(label = ifelse(pos.lfc >= 2.2, sgrna_gene, "")),
                  position = position_jitter(seed = 31),
                  size=3, box.padding = 1)+
        scale_color_manual(values =c("orange", "#006400", "grey")) +
        #scale_color_manual(values =c("blue","orange", "red", "grey")) +
        scale_size(range = c(0.01, 10)) +
        theme_bw()+theme(panel.grid =element_blank())+labs(title="Donor3, Count ≥ 30")
#ggsave(plot=p_donor3, filename="/home/chengennong/ylp/aav/fig_CapIV/donor_d3.png", dpi = 300, width=9,height=5)
ggsave(plot=p_bubble, filename="/home/chengennong/ylp/aav/screen/screen_out/batch4/plot/1a.png", dpi = 300, width=6,height=5)
#ggsave(plot=p_donor3_text, filename="/home/chengennong/ylp/aav/fig_CapIV/ye/donor_d3_text.png", dpi = 300, width=9,height=5)

# venn
a_marker_sgrna_gene_tests_table %>% arrange(-pos.lfc) %>% head(200) %>% pull(sgrna_gene)

vec_aav2 <- c("KIAA0319L", "VPS29", "GPR108", "TM9SF2", "B3GAT3", "VPS52", "VPS54", "SLC35B2", "ATP2C1", "RNF121", "EXT1", "JTB", "RGP1", "VPS53", "SWIP", "VPS51", "B4GALT7", "XYLT2", "ARPC2", "NDST1", "OSBPL11", "B3GALT6", "RIC1", "C16orf62", "COG8", "ATP6V0A2", "NAA38", "OSBPL9", "EXT2", "ARPC4", "VPS35", "COMMD3", "HTT", "CCDC22", "DMXL1", "TRAPPC13", "RAB6A", "ACP2", "GPR107", "MALAT1", "SPG8", "COG7", "RAB7A", "RABIF", "PBRM1", "COMMD8", "FLJ37453", "ZBTB32", "ACTR2", "TSPAN4", "SPOP", "KAZALD1", "OGFOD2", "RHBDD3", "DKFZP686I15217", "ARID2", "TICAM1", "HOXA10", "NCRNA00167", "TAF12", "PAX2", "C5orf39", "C3orf67", "SEC14L2", "COG4", "HGS", "WRB", "SIX3", "TCF19", "CCDC93", "PLEKHG2", "RGL2", "OSM", "MDK", "OMG", "STOX2", "CCDC53", "COG1", "SLC39A7", "TP53I13", "GTF2I", "KDM2B", "DNAJC30", "C14orf176", "SCFD1", "MKX", "C14orf4", "LOC100288730", "TSSC1", "RAB2A", "CIZ1", "INCA1", "SERPINB7", "FBXW5", "GRRP1", "DOPEY2", "NOL3", "SPHK1", "C22orf34", "SND1", "ZNF436", "COG2", "DLX4", "P2RY2", "XPO1", "SP1", "C3orf71", "FAM69B", "RDH11", "TCEA2", "STRN4", "FLAD1", "MTHFSD", "NFXL1", "PDIK1L", "NFIX", "COG3", "CRTAM", "LTBP1", "SDC1", "YAF2", "ZMYM3", "F12", "GLG1", "PPIL3", "WDR20", "BRD2", "C7orf50", "LRRC14", "CLEC16A", "C9orf95", "USP47", "ZNF628", "ZSWIM3", "SRP72", "UNK", "DNAJC27", "RFXANK", "STC2", "TMUB1", "ATF4", "AZI2", "C11orf71", "CENPB", "CLDN22", "FAM187B", "HIVEP1", "ILKAP", "PRKX", "RAB1A", "SLC35D1", "TMEM208", "VPS39", "ZNF268", "NYX", "ARRDC4", "ATP6V1B1", "ATP6V1C1", "BZRAP1", "C7orf23", "CLN8", "CROCC", "CWC25", "CXXC4", "DNAJB4", "FAAH2", "FAM127A", "GALNT12", "GJD3", "GPRASP2", "HAND2", "HSF2", "ING5", "IQCE", "JMJD8", "LYVE1", "NAGA", "NCRNA00086", "NHLRC4", "P2RX7", "PIGM", "PLA2G3", "SNORD26", "TAS2R39", "THSD1P1", "TREX1", "VPS37D", "ZDHHC6", "ZNF37B", "FASTKD3", "FLJ16779", "GATAD1", "HCFC1", "LEPREL2", "LOC400657", "NAA35", "SCAND1", "SLC17A9")
intersect(vec_aav2, a_marker_sgrna_gene_tests_table %>% arrange(-pos.lfc) %>% head(1000) %>% pull(sgrna_gene)) %>% sort