#env cnmf
######################
ada_1=GGCTGTCTCTTATACACATCTCCGA
ada_2=GCTCTGTCTCTTATACACATCTGAC
anchor_1_3end=TTGCTGTTTAGCCGTGGGTCTCCAGCTGGCATGTCT
anchor_2_3end=GGACTGATTTTGAGTTCTGTTCAGGTAA
anchor_1_5end=agcTTACCTGAACAGAACTCAAAATCAGTCC
anchor_2_5end=cCAGACATGCCAGCTGGAGACCCACGGCTAAACAGCAA
batch_dir=/mnt/cgn/data/aav/aav_batch3_2023_1218 #/mnt/cgn/data/aav/aav_batch1_2023_1127 /mnt/cgn/data/aav/aav_batch3_2023_1218 /mnt/cgn/data/aav/aav_batch2_2023_1206
for sample in LH007 LH008 LH009 LH010 LH011 LH012 LH013 #LH001 LH002 LH003 LH004 LH005 LH006 
do
    cd ${batch_dir}/Rawdata/${sample}; echo ${sample}
    #cut adaptor
    cutadapt -a ${ada_1} -A ${ada_2} -j 24 -e 0.1 \
        -o ${sample}_R1.cut.fq.gz  -p ${sample}_R2.cut.fq.gz  ${sample}_R1.fq.gz  ${sample}_R2.fq.gz > cut_adaptor.log
    #quality control: q>25
    sickle pe -q 25 -t sanger -g -f ${sample}_R1.cut.fq.gz -r ${sample}_R2.cut.fq.gz \
        -o ${sample}_R1.sickle.fq.gz -p ${sample}_R2.sickle.fq.gz -s ${sample}_single_trim.sickle.fq.gz > sickle.log
    #trim 3' arm
    cutadapt -a ${anchor_1_3end} -A ${anchor_2_3end} -j 24 -e 0.1 \
        -o ${sample}_R1.anchor3.fq.gz  -p ${sample}_R2.anchor3.fq.gz  ${sample}_R1.sickle.fq.gz  ${sample}_R2.sickle.fq.gz > cut_anchor3.log
    #trim 5' arm
    cutadapt -g ${anchor_1_5end} -G ${anchor_2_5end} -j 24 -e 0.1 \
        -o ${sample}_R1.anchor53.fq.gz  -p ${sample}_R2.anchor53.fq.gz  ${sample}_R1.anchor3.fq.gz  ${sample}_R2.anchor3.fq.gz > cut_anchor5.log
    #filter reads with N base > 0 and keep only 21 nt
    cutadapt -j 24 -m 21 -M 21 --max-n 0 \
        -o ${sample}_R1.7mers.fq.gz  -p ${sample}_R2.7mers.fq.gz  ${sample}_R1.anchor53.fq.gz  ${sample}_R2.anchor53.fq.gz > cut_7mers.log
    seqkit stat ${sample}*.fq.gz > ${sample}.stat
    zcat ${sample}_R1.7mers.fq.gz|awk '(NR-2)%4 == 0'|sort|groupBy -g 1 -c 1 -o count| sort -nrk2 | \
        awk -v sample=$sample 'BEGIN {print "Library\tSequence\tCount"} {print sample"\t"$0}' > ${sample}_R1.7mers.nt.stat
    zcat ${sample}_R1.7mers.fq.gz|seqkit translate|awk '(NR-2)%2 == 0'|sort|groupBy -g 1 -c 1 -o count| sort -nrk2 | \
        awk -v sample=$sample 'BEGIN {print "Library\tSequence\tCount"} {print sample"\t"$0}' > ${sample}_R1.7mers.aa.stat
    cd ${batch_dir}/Rawdata
done


