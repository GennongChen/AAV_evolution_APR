

cd /home/chengennong/ylp/aav/rnaseq
snakemake --snakefile /home/chengennong/ylp/aav/rnaseq/run_bulk_pre.smk \
    --cores 78 --rerun-incomplete  -n
cd /mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata
cut -f1,7- counts.txt | awk 'NR>1' > count1.txt