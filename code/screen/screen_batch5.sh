
#cda genome
seqkit=/home/chengennong/tools/mambaforge/envs/cnmf/bin/seqkit
cutadapt=/home/chengennong/tools/mambaforge/envs/cnmf/bin/cutadapt
fastp=/home/chengennong/tools/mambaforge/envs/cnmf/bin/fastp
bowtie2=/home/chengennong/tools/bowtie2/bowtie2
samtools=/home/zhanghaosheng/mydisk/miniconda3/envs/rserver/bin/samtools
bedtools=/usr/bin/bedtools

raw_dir=/mnt/cgn/data/aav/screen_${batch}_2025_0507/Rawdata
batch=batch5
#for i in `seq 10 24`
#do mkdir -p ${raw_dir}/Evo${i}
#done

#qc
cd ${raw_dir}
for i in `ls`
do echo $i && cd $i
    R1_fq=${i}_R1.fq.gz
    R2_fq=${i}_R2.fq.gz
    $fastp \
    -i $R1_fq \
    -I $R2_fq \
    -o ${i}_R1.fastq.gz \
    -O ${i}_R2.fastq.gz \
    --qualified_quality_phred 20 \
    --thread 16 > fastp.log
    cd ${raw_dir}
done

#demulti & cutadapt
screening_meta_table=/home/chengennong/ylp/aav/screen/screen_${batch}_sample_info.txt
for i in `awk 'NR>1' $screening_meta_table |cut -f1|sort -u`
do echo $i && cd ${raw_dir}/$i
    #i5_index=`awk 'NR>1' $screening_meta_table | grep $i |cut -f7|sort -u`
    #input_name=`awk 'NR>1' $screening_meta_table | grep $i |cut -f9|sort -u`
    ada5=`awk 'NR>1' $screening_meta_table |cut -f10|sort -u`
    ada3=`awk 'NR>1' $screening_meta_table |cut -f11|sort -u`
    #$cutadapt -j 20 -g ^${i5_index} -o ${i}_R11.fastq.gz ${raw_dir}/$input_name/R1_qc.fastq.gz
    $cutadapt -j 20 -a $ada3 --error-rate 0.2 -o ${i}cut3_R1.fastq.gz ${i}_R1.fastq.gz > cut3.log
    $cutadapt -j 20 -g $ada5 --error-rate 0.2 -o ${i}cut35_R1.fastq.gz ${i}cut3_R1.fastq.gz > cut5.log
done
ll * 

#count
cat /home/chengennong/ylp/aav/screen/sgRNA_library/Broadgpp.txt | cut -f7,2 | \
    sed 's/Non-Targeting Control/NTC/g' | \
    awk 'NR>1 {print $1"_"$2"\t"$2"\t"$1}' > /home/chengennong/ylp/aav/screen/sgRNA_library/CRISPRko.Broadgpp.txt.library

library_dir=/home/chengennong/ylp/aav/screen/sgRNA_library
out_dir=/home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count
mageck=/home/chengennong/tools/mambaforge/envs/genome/bin/mageck
mkdir -p $out_dir && cd ${raw_dir}
screening_meta_table=/home/chengennong/ylp/aav/screen/screen_${batch}_sample_info.txt


for Library_i in `awk 'NR>1' $screening_meta_table |cut -f3|sort -u`
do
    grep "$Library_i" $screening_meta_table
    library=${library_dir}/`grep "$Library_i" $screening_meta_table|cut -f2,3|sed 's/\t/./'|sort -u`.txt.library
    out_prefix=${out_dir}/`grep "$Library_i" $screening_meta_table|cut -f2,3|sed 's/\t/./'|sort -u`
    sample_label=`cat $screening_meta_table|sort -k5,6 |grep "$Library_i"|cut -f9|sed 's/\t/_/g'|xargs echo|sed 's/ /,/g'`
    fastq1=`cat $screening_meta_table|sort -k5,6 |grep "$Library_i"|cut -f1|xargs -I {} echo ${raw_dir}/{}/{}cut35_R1.fastq.gz`
    #fastq2=`grep "$Library_i" $screening_meta_table|cut -f1|awk '{print  "${raw_dir}/" $1 "/" $1 "_R2.fastq.gz"}'|xargs echo`
    #echo $Library_i $library $out_prefix $sample_label $fastq1
    $mageck count \
    -l $library \
    -n $out_prefix \
    --sample-label $sample_label \
    --fastq $fastq1 \
    --sgrna-len 20 \
    --norm-method none #--trim-5 22,23,24,25,26,28,29,30
done
head /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.count.txt
grep -Ew "CD7|KIAA0319L" /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.count.txt #|sort -k 12,12n|cut -f1,12-

#2 Next, raw read counts across both library sets were normalized to the total read count in each sample, and each of the matching samples across two sets were merged to generate a single normalized read count table.
/home/chengennong/tools/mambaforge/envs/genome/bin/Rscript \
    /home/chengennong/ylp/aav/code/screen_count_nor.R \
    -i /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.count.txt \
    -o /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.normalized_count.txt \
    -n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.no_target.txt
#3 test #cda genome
$mageck test \
    -k /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.normalized_count.txt \
    -t Cell-1-P4,Cell-2-P4,Cell-3-P4 \
    -c Cell-1-P7,Cell-2-P7,Cell-3-P7 \
    --norm-method none \
    --control-sgrna /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.no_target.txt \
    -n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7 --gene-lfc-method median --paired
#$mageck test \
#    -k /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.normalized_count.txt \
#    -t Cell-1-P4,Cell-2-P4,Cell-3-P4 \
#    -c Cell-1-P6,Cell-2-P6,Cell-3-P6 \
#    --norm-method none \
#    --control-sgrna /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.no_target.txt \
#    -n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P6 --gene-lfc-method median --paired
$mageck test \
    -k /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.normalized_count.txt \
    -t Cell-1-P4,Cell-2-P4,Cell-3-P4 \
    -c Cell-1-baseline,Cell-2-baseline,Cell-3-baseline \
    --norm-method none \
    --control-sgrna /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.no_target.txt \
    -n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_baseline --gene-lfc-method median --paired

#neg
#pos
#sort -k 12,12n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_unsorted.gene_summary.txt | head -n15
sort -k 12,12n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7.gene_summary.txt | head -n15 |cut -f1,12-
sort -k 14,14rn /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7.gene_summary.txt | head -n15|cut -f1,12-|grep -v "e-"

sort -k 12,12n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P6.gene_summary.txt | head -n15 |cut -f1,12-
sort -k 14,14rn /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P6.gene_summary.txt | head -n15|cut -f1,12-|grep -v "e-"

sort -k 12,12n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_baseline.gene_summary.txt | head -n15 |cut -f1,12-
sort -k 14,14rn /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_baseline.gene_summary.txt | head -n15|cut -f1,12-|grep -v "e-"

grep -Ew "CD7|TM9SF2|VPS52|KIAA0319L|NR4A1|KMT2D|RUNX3" /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7.gene_summary.txt |sort -k 12,12n|cut -f1,12-

sort -k 12,12n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7.gene_summary.txt|awk '$14 > 1' | wc -l
#plot
$mageck plot -k /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.normalized_count.txt \
    -g /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7.gene_summary.txt \
    --norm-method none \
    --control-sgrna /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/CRISPRko.Broadgpp.no_target.txt \
    -n /home/chengennong/ylp/aav/screen/screen_out/${batch}/mageck_count/Broadgpp_P4_P7 \
    --genes CD7,TM9SF2,VPS52,KIAA0319L,NR4A1,KMT2D,RUNX3,GPR108,NTC

mkdir /home/chengennong/ylp/aav/screen/screen_out/batch5/plot
for sc in `ls /home/chengennong/ylp/aav/screen/screen_out/fig_script`
do #/home/chengennong/tools/mambaforge/envs/seurat4/bin/Rscript \
    echo /home/chengennong/ylp/aav/screen/screen_out/fig_script/$sc
done
/home/chengennong/tools/mambaforge/envs/seurat4/bin/Rscript /home/chengennong/ylp/aav/screen/screen_out/batch5/fig_script/fig1c-d.R
#4 Gene hits were classified as having a median absolute log2-fold change >0.5 and a false discovery rate (FDR) <0.05. For supplemental CD4+ screens (fig. S9), reads were aligned to the full Calabrese A and B library in a single reference file.
#statistic & draw result
#for i in `ls /home/chengennong/code-manual/vscode/ylp/Science2022/script_dir/fig_script_dir/* |grep -v "#" |grep fig[12]`
#do /home/chengennong/anaconda3/envs/R/bin/Rscript $i
#done


