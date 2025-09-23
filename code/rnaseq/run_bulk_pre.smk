
# Snakefile
import pandas as pd

# 读取样本ID列表
samples_df = pd.read_csv("/home/chengennong/ylp/aav/rnaseq/metatable.txt", sep="\t")  # 假设分隔符是制表符，也可以用逗号或者空格
SAMPLES = samples_df['Sample'].tolist()

FASTP = "/home/chengennong/tools/mambaforge/envs/cnmf/bin/fastp"
STAR = "/home/chengennong/tools/mambaforge/envs/cnmf/bin/STAR"
FEATURECOUNTS = "/home/chengennong/tools/mambaforge/envs/cnmf/bin/featureCounts"


# 设定输入样本
#SAMPLES, = glob_wildcards("data/{sample}_R1.fq.gz")

# 目标文件，最终输出的基因计数结果
TARGETS = ["counts.txt"]

# work_dir
workdir: "/mnt/cgn/data/aav/rnaseq_batch1_2025_0107/Rawdata"

# 运行所有规则
rule all:
    input:
        TARGETS

# 规则: 质控
rule fastp:
    input:
        r1="{sample}/{sample}_R1.fq.gz",
        r2="{sample}/{sample}_R2.fq.gz"
    output:
        r1_clean="{sample}/{sample}_R1.cut.fq.gz",
        r2_clean="{sample}/{sample}_R2.cut.fq.gz",
        html="{sample}/{sample}_fastp.html",
        json="{sample}/{sample}_fastp.json"
    threads: 6
    params:
        n_base_limit=10,
        qualified_quality_phred=15,
        unqualified_percent_limit=40
    shell:
        """
        {FASTP} --thread {threads} \
            --n_base_limit {params.n_base_limit} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --unqualified_percent_limit {params.unqualified_percent_limit} \
            -i {input.r1} -I {input.r2} \
            -o {output.r1_clean} -O {output.r2_clean} \
            --html {output.html} --json {output.json}
        """

# 规则: 使用STAR进行比对
rule star_alignment:
    input:
        r1="{sample}/{sample}_R1.cut.fq.gz",
        r2="{sample}/{sample}_R2.cut.fq.gz"
    output:
        bam="{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    params:
        genome_dir="/mnt/cgn/data/ref/human/refdata-gex-GRCh38-2020-A/star/",
        read_files_command="zcat",
        transcriptome_sam="TranscriptomeSAM",
        mismatches=10
    threads: 6
    shell:
        """
        {STAR} --runThreadN {threads} --runMode alignReads \
            --outFilterMismatchNmax {params.mismatches} \
            --readFilesCommand {params.read_files_command} \
            --quantMode {params.transcriptome_sam} \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped None \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {wildcards.sample}/{wildcards.sample}.
        """

# 规则: 使用featureCounts进行定量
rule featurecounts:
    input:
        bam=expand("{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
    output:
        counts="counts.txt"
    params:
        gtf="/mnt/cgn/data/ref/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        threads=20
    shell:
        """
        {FEATURECOUNTS} -t exon -g gene_id -T {params.threads} \
            -p --countReadPairs -O -g gene_name \
            -a {params.gtf} -o {output.counts} {input.bam}
        """
