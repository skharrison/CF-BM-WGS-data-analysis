
# Burkholderia WGS Analysis 

 - Repository containing all pieces and code involved in analyzing whole genome sequences of 19 Burkholderia multivorans isolates over 3 year time frame from sputum samples of a patient with cystic fibrosis.  


### Data Quality Control
- Raw reads trimmed with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) using the following parameters:
``` 
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:80
```
- Visualized reads using [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) manually and trimmed more as needed

### Reference Mapping and Variant Calling
- Trimmed reads mapped using [BWA](http://bio-bwa.sourceforge.net/) to generate SAM file. Duplicates marked and sorted by [GATK4](https://gatk.broadinstitute.org/hc/en-us) and [Samtools](http://www.htslib.org/).
- Script snippet to accomplish this:
```
bwa mem -M -R '@RG\tID:HSQ-7001360\tSM:AS${sample}\tLB:library1\tPL:ILLUMINA' -M -t 16 GCF_003019965.1_ASM301996v1_genomic.fna ~/trimmed_wgs3/AS${sample}_R1_001_paired.fastq.gz ~/trimmed_wgs3/AS${sample}_R2_001_paired.fastq.gz > AS${sample}.sam
samtools fixmate -O bam AS${sample}.sam AS${sample}.bam
samtools sort -O bam -o AS${sample}.sorted.bam AS${sample}.sam
gatk MarkDuplicatesSpark -I AS${sample}.sorted.bam -O AS${sample}_marked_duplicates.bam
samtools index AS${sample}_marked_duplicates.bam
```
- Samtools used to call mpileup files:
```add mpileup here
```
- Code written to call variants from mpileups:
```
add code here can make cell python color coding 
```

- Wasing thinking after this could then add link to send to a jupyter notebook that contains all code for everything like allele frequency cutoff generation
- Pandas table and manipulation
- R used to generate final variants kept based on standard deviation 

**VCF annotation/analysis**
- code for vcf annotation (perl/python?), KEGG/KOG protein prediction
- script to remove common variants to all (bcftool commands)
- scripts and code for extraction of desired info
- scripts for plotting and visualization (notebooks?)

**Phylogenetic Analysis**
- SNP alignment (software, scripts, code) .. possibly just whole genome alignment??
- raxml script (script, options chosen)
- visualization (software, options chosen, outgroup/rooting choices)

**Related information**
- Background info on CF/Burkholderia/bacterial evolution in vivo
- bacterial strain information/tables 
- methods for sequencing/isolation/culturing


