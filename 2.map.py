#gatk HaplotypeCaller -R ~/database/reference/GCF_000001405.40_GRCh38.p14_genomic.fasta $(ls *sorted.bam | awk '{printf "-I %s ", $0}') -O out_test.vcf

import os

from glob import glob

reference = '/home/cat/database/reference/GCF_000001405.40_GRCh38.p14_genomic.fasta'
# Map the paired-end reads
for file in glob("./data/*_trimmed_R1*"):
    # Extract the base of the filename
    base = file.split("_trimmed_R1")[0]
    
    # Map the reads to the reference genome
    rgid = base.split('/')[-1].split('_')[0]
    com = fr"bwa mem -t 12 {reference} {base}_trimmed_R1.fastq.gz {base}_trimmed_R2.fastq.gz -R '@RG\tID:{rgid}\tSM:{rgid}\tLB:Library1\tPL:illumina\tPU:PlatformUnit1' |samtools sort -o {base}.sorted.bam"
    os.system(com)

    # index BAM file
    os.system(f"samtools index {base}.sorted.bam")

