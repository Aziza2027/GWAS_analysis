import os

from glob import glob

reference = '/home/cat/database/reference/GCF_000001405.40_GRCh38.p14_genomic.fasta'
# Map the paired-end reads
for file in glob("./data/*_trimmed_R1*"):
    # Extract the base of the filename
    base = file.split("_trimmed_R1")[0]
    
    # Map the reads to the reference genome
    os.system(f"bwa mem {reference} {base}_trimmed_R1.fastq.gz {base}_trimmed_R2.fastq.gz |samtools sort -o {base}.sorted.bam")

    # index BAM file
    os.system(f"samtools index {base}.sorted.bam")

