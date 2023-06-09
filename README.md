# GWAS_analysis

![IM](im.png)


[Image source](https://www.researchgate.net/figure/Schematic-of-a-bioinformatics-pipeline-Examples-of-the-most-commonly-used-publicly_fig3_250923605)

## *What is Bioinformatics pipeline?*

*Bioinformatics pipelines are used to process large amounts of biological data, such as sequencing data, and to extract meaningful insights from it. The steps involved in a bioinformatics pipeline can vary depending on the specific analysis being performed, but here are some general steps that are typically involved:*

 - *Step 1: Read Preprocessing During this step, the raw sequencing data
   is processed to prepare it for analysis. Tasks such as trimming
   adapter sequences, filtering out low-quality reads, and removing any
   contaminants are performed.*
   
 - *Step 2: Alignment In this step, the preprocessed reads are aligned to
   a reference genome or transcriptome. Alignment algorithms are used to
   identify the positions in the reference where the reads originated
   from.*
   
  - *Step 3: Variant Calling After alignment, the next step is variant
   calling. This involves identifying genetic variations, such as single
   nucleotide polymorphisms (SNPs), insertions, or deletions, by
   comparing the aligned reads to the reference sequence.*
   
   - *Step 4: Annotation Once the variants have been identified, they are
   annotated to determine their potential functional effects. Annotation
   tools are used to provide information about the location of the
   variants in the genome, their impact on genes, and potential disease
   associations.*
   
   - *Step 5: Fill Missing Genotypes In this final step, any missing
   genotypes in the variant data are inferred or imputed. This can be
   done by leveraging information from neighboring variants or using
   statistical methods to estimate the missing genotypes based on the
   available data.*

*It's important to note that these steps are commonly followed in genomic data analysis pipelines, but the specific tools and techniques used may vary depending on the analysis goals, the type of data, and the available resources.*

Text was generated using [ChatGPT](https://chat.openai.com/)

Requirements:

 - python 3.8-3.10 
 Scikit-learn, numpy, pandas 
 - Fastp 
 - Samtools 
 - BWA
 - Bcftools
 - tabix
 - reference genome + dbSNP
## Step 1(Read Preprocessing)

Put all of your *fastq.gz files to data folder and run:
-     python 1.trim.py

It removes low-quality reads and trims adapter sequences. It also trims reads that are longer than the desired length and filters out reads that are too short and generates reports.


## Step 2(Alignment)

Make sure your data folder contains *trimmed* files and run:
-     python 2.map.py

This command aligns reads to a reference genome and generates .bam files(+.bai files).
*Do not forget to download, index [reference genome](https://www.ncbi.nlm.nih.gov/genome/guide/human/) and change path.*

## Step 3(Variant Calling)

You may have to change permission before running `3.var_call` and make sure configurations meet your requirements. `AN` means Total number of alleles in called genotypes.
Run following command:

-     ./3.var_call

Output is saved in `variants` folder(separate files for `SNPs` and `INDELs`)

## Step 4(Annotation)

Download [dbSNP](https://www.ncbi.nlm.nih.gov/genome/guide/human/) and index(maybe there is easier way, but that's how I did).
Run following command:

-     ./4.annotate_var

Output is saved in `variants` folder(annotates only SNPs not INDELs, you have to do some modifications to annotate INDELs).

## Step 5(Fill Missing Genotypes)

Run:
-     python 5.fill_missing.py

It fills missing values. Output files are saved in `variants` folder. `original.csv` original file with variants, `filled.csv` file with imputed genotypes.


### Links


[Website](https://aziza-diabetes.netlify.app) 


[Video Tutorial](https://youtu.be/RSbNVv47EHM) 


[GWAS results(INDELS)](https://drive.google.com/file/d/1eohdtdyZrF-Xog0zflciWzuxy2kYq7Cm/view?usp=sharing) 

[GWAS results(SNPs)](https://docs.google.com/spreadsheets/d/1OIdxdldcfwFWI-uHqFqgwsSimjU7js34/edit?usp=share_link&ouid=112267455847643950553&rtpof=true&sd=true) 

[GWAS results(Hardy-Weinberg Equilibrium)](https://drive.google.com/file/d/1ephEzd4ung1l5aFKupGyZW5X5AIgxI8s/view?usp=share_link) 

[GWAS results(P-value)](https://drive.google.com/file/d/12k5GoLTlckXjSBJef40jBR2oq2phNwdY/view?usp=sharing) 

[GWAS results(OR,CI)](https://drive.google.com/file/d/1PKaDUQjyc6-NDEtB8oMxD-kHix50X1PS/view?usp=sharing) for SNPs with P-value less than 0.05(codominant model).

[GWAS results(P-value)](https://docs.google.com/spreadsheets/d/1lwYOIycMqoL9Km2Sd4-qq185JmZd25Iobf3VUXnWejs/edit?usp=drive_link) without HWE test

[Excluded sampes + SNPs(2+alleles)](https://docs.google.com/document/d/1Q7PxEzHJMG7UideHzl2PTs6Em3BkV62mrLiwfAMiDR4/edit?usp=sharing) 




