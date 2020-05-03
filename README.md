
# Burkholderia WGS Analysis 


* Repository containing all pieces and code involved in analyzing whole genome sequence short read data of 17 Burkholderia multivorans isolates over 3 year time frame from sputum samples of a patient with cystic fibrosis.
* Highly recommend use of [Anaconda](https://www.anaconda.com/products/individual) to install software components.  
* In order to carry out SNP analysis, interactive python3 [Jupyter notebook](https://jupyter.org/) utilized, also easily installed via Anaconda. 

**NOTE:** ALL samples analyzed were sequenced with 151bp Illumina paired end reads  

### Data Quality Control
- Raw reads trimmed with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) using the following parameters:
``` 
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:100
```
> - ILLUMINACLIP: Contained TruSeq3 primers.  
> - LEADING/TRAILING: Removing poor qulaity beginning and ends. 
> - SLIDINGWINDOW: If any window of 4 reads quality drops below 15 will cut sequence. 
> - MINLEN: Chose a conservative 100 to ensure only analyzing high quality read data.  

- Visualized reads using [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) manually and trimmed more as needed.

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
- [Samtools](http://www.htslib.org/) used to call mpileup files:
```add mpileup command here
list_of_samples="218 219 222 223 224 225 226 227 228 229 230 231 232 233 236 237 240"
for sample in $list_of_samples; do
samtools mpileup -f GCF_003019965.1_ASM301996v1_genomic.fna -B -R -aa AS${sample}.sorted.bam -o AS${sample}.mpileup
done
```
[Script: parse_mpileups.py](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/scripts_notebooks/parse_mpileups.py) grabs single nucleotide variants positions from mpileup and generates a table for each strain containing alternate allele and allele frequency (# alt alleles at position)/(total reads at position). 

[Jupyter Notebook: parse_mut_calls.ipynb](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/scripts_notebooks/parse_mut_calls.ipynb)

**Contains code to:**
- Remove frequencies that fall below 0.03 in each table (remove variants due to sequence error)
- Merge all strain tables to one large dataframe
- Keep sites that rise or fall in frequency (Removing variant positions common to all strains with >.95 freq)

Table Example:

![merged_table](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/table_image.png)

### SNP Position Analysis 
-------
[Jupyter Notebook: SNP_Analysis.ipynb](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/scripts_notebooks/SNP_Analysis.ipynb) 
  
**Contains code to:**
- Plot distribution of filtered SNP locations across all chromosomes
- Extract all genes located in SNP rich locations from GFF file (from same reference that mapped reads to)
- Extract all genes with at least one SNP and annotate with protein product when CDS, also count number of SNPs per gene using gene start/stop locations and variant position numbers

- Table created displaying top 15 protein coding genes with highest numbers of SNPs:

![Table image](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/table_image.png)


### Genome Assembly
- Used [Spades](http://home.cc.umanitoba.ca/~psgendb/doc/spades/manual.html) on all previously trimmed paired end data. Now have files such as, AS{sampleID}_R1.fastq and AS{sampleID}_R2.fastq. Then used bash script with loop to accomplish all files. Bash script snippet:

**INPUT FILE FORMAT**: sampleID_R1.fastq, sampleID_R2.fastq .   
**INPUUT FILE EXAMPLE**: AS218_R1.fastq, AS218_R2.fastq .   

```
for file1 in /path/to/trimmed/data/directory/*R1*fastq
do 
file2=${file1/R1/R2}
out=${file1%%_R1.fastq}_output
spades.py -k 21,33,55,77 -c --only-assembler -1 $file1 -2 $file2 -o $out
done
```
> -c minimizes mismatches.   
> -k 21,33,5,77 because recommended to use for short paired end Illumina reads.    
> -1 forward read file -2 reverse read file.    
> -o output file that specified name by replacing the _R1.fastq portion of input file with _output. Allowing for each sample to then be in its own directory.    

- Annotated genomes using Prokka (add command used below)

```
prokka command here
```

### Genome Comparative Analysis and Rearrangement Analysis
- Using [Mauve](http://darlinglab.org/mauve/download.html) to generate whole genome alignments to look at recombinant regions.  

**IMPORTANT NOTE:** Mauve only accepts genbank input as .gbk extension, also accepts fasta format so can use reference fasta and contigs obtained from assembly but if do want use genbank format can easily change extension using bash loop:

```
for file in *.gbff; do
mv -- "$file" "${file%.gbff}.gbk"
done
```
- [ ] Order contigs using progressive mauve
- [ ] Align all genomes to one another 
- [ ] Potentially use [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML)

### Phylogenetic Analysis
STEPS PLAN TO DO:
- [ ] Need to do recombinant analysis first as those regions can interfere with phylogeny result 
- [ ] Obtain high quality variants for each strain 
- [ ] Extract all genes containing variants from reference using Bipython
- [ ] Write python code to modify reference gene sequences for each sample to contain ALT positions
- [ ] Concat all genes together to create consensus fasta sequence for each sample (header being sample name)
- [ ] Concat all samples together to create one multifastA 
- [ ] Use [MUSCLE](https://biopython.org/DIST/docs/api/Bio.Align.Applications._Muscle.MuscleCommandline-class.html) or [MAFFT](https://biopython.org/DIST/docs/api/Bio.Align.Applications._Mafft.MafftCommandline-class.html) python wrapper to create sequence alignment 
- [ ] Input sequence alignment into RAXML using following command line:

```
raxmlHPC -f a -x 43734 -p 89493 -# 100 -s snp_alignment.fasta -n phylo_result -m GTRGAMMA 
```
- [ ] Visualize tree using [FigTree](http://evomics.org/resources/software/molecular-evolution-software/figtree/) 

### Automating steps??? (get to this hopefully eventually)
- Snakemake pipeline one day 

