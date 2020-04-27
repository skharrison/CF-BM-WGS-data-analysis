
# Burkholderia WGS Analysis 

 - Repository containing all pieces and code involved in analyzing whole genome sequences of 19 Burkholderia multivorans isolates over 3 year time frame from sputum samples of a patient with cystic fibrosis.
 - potentially add beginning step to install using conda at each step for each piece of software (all is definetely on conda) just to make super easy to implement if make code snippets more transferrable.
 - also somewhere would like to add a bash script template that can easily be modified for this project that allows everything to work when submitting jobs to a queing system. 



### Data Quality Control
- add in script with code snippet to accomblish for x # of samples 
- Raw reads trimmed with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) using the following parameters:
``` 
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:100
```
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
- Samtools used to call mpileup files:
```add mpileup command here
```

### Genome Variation/Mutation Analysis
- Python code either link to file or below to get sequence error cutoff 
- Python code to get extract variants and generate freq tables etc..  





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
> -k 21,33,5,77 because recomennded to use for short paired end Illumina reads.    
> -1 forward read file -2 reverse read file.    
> -o output file that specified name by replacing the _R1.fastq portion of input file with _output. Allowing for each sample to then be in its own diretory.    

- Annotated genomes using Prokka (add command used below)

```
prokka command here
```

### Genome Comparative Analysis and Rearrangement Analysis
- Using [Mauve](http://darlinglab.org/mauve/download.html) to conduct comparisons against complete Burkholderia _multivorans_ on refseq [NCBI](https://www.ncbi.nlm.nih.gov/assembly). Downloaded all complete genomes in Genbank gbff format (there is 15).  

**IMPORTANT NOTE:** Mauve only accepts genbank input as .gbk extension so a helpful bash loop to rename all file extensions in a directory

```
for file in *.gbff; do
mv -- "$file" "${file%.gbff}.gbk"
done
```
- add code for generating alignments via progressive mauve
- add all steps for using ClonalOrigin to generate rearrangment maps to be analyzed

### Phame Phylogenetic Analysis (maybe..)


### Automating steps??? (get to this hopefully eventually)
- Snakemake pipeline one day 

### Helpful Bioinformatics along the way... 
- Add cool/helpful places that have found helpful while doing bacterial genome analysis 
