**By: Sarah Harrison and John Navarro** 


# Burkholderia WGS Analysis 


* Repository containing all pieces and code involved in analyzing whole genome sequence short read data of 17 Burkholderia multivorans isolates over 3 year time frame from sputum samples of a patient with cystic fibrosis.
* In order to carry out SNP analysis, interactive python3 [Jupyter notebook](https://jupyter.org/) utilized, installed via Anaconda. 

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
samtools index AS${sample}.sorted.bam
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

![merged_table](https://github.com/skharrison/CF-BM-WGS-data-analysis/blob/master/merged_table.png)

### SNP Position Analysis


### Structural Variant Analysis 

**Programs used:**
- [Pindel](https://github.com/genome/pindel) Helpful documentation [here](http://gmt.genome.wustl.edu/packages/pindel/user-manual.html)
- [Breakdancer](https://github.com/genome/breakdancer) Helpful documentation [here](https://gmt.genome.wustl.edu/packages/breakdancer/documentation.html) 
- [DELLY](https://github.com/dellytools/delly)
- [Manta](https://github.com/Illumina/manta)
- [LUMPY](https://github.com/arq5x/lumpy-sv)
- [Gridss](https://github.com/PapenfussLab/gridss)

pindel commands:
```
## first have to generate config files for each sample, format of file shown below 
## created each samples pindel configs by using samtools stats bam and using average insert size (saved each config as sample{sampleNum}.txt)
path/to/bam avginsertSize sampleName

## also created ploidy file that contained each chrom name and the ploidy (since bacterial samples ploidy 1 next to each chrom name)

# run pindel 
pindel -f ref.fna -i sample_config.txt -o p1 -c ALL -T 8 -Y -a 4 -m 4 -Y ploidy.txt 

# make vcf of output 
/usr/local/pindel/pindel2vcf -P p1 -r ref.fna -R SALREF -d 20200823 -v sample_pindel.vcf
```
breakdancer commands:
```
# config creation:
perl /usr/local/breakdancer/perl/bam2cfg.pl -g -h ~/sample.bam > sample.cfg
# run breakdancer: 
breakdancer-max -r 4 sample.cfg > sample_breakdancer.vcf
```

delly commands:
```
# run delly on each bam 
delly call -g ref.fna -o ample.bcf sample.bam

# after all samples finished merge together
delly merge -o all_sites.bcf *.bcf

# genotype each sample from merged 
delly call -g ref.fna -v all_sites.bcf -o /sample$_delly.bcf sample.bam

# convert bcf to vcf 
bcftools view sample_delly.bcf > sample_delly.vcf
```

manta commands:
```
# NOTE: made sample directories in manta directory prior to running by:
for num in {1..17} ; do 
mkdir sample${num}
done

#set everything up
configManta.py --bam sample.bam --referenceFasta ref.fna --runDir /sample

#run workflow inside sample directory
/sample/runWorkflow.py

```
lumpy commands:
```
##generates insert size statistics and generates needed .histo file 
samtools view sample.bam \
    | tail -n+100000 \
    | /usr/local/lumpy-sv/scripts/pairend_distro.py \
    -r 151 \
    -X 4 \
    -N 10000 \
    -o sample.lib1.histo
——> this script gives mean/std output that looks like:
Removed 6 outliers with isize >= 1328
mean:311.452552553	stdev:125.998418651

— Used the mean and stdev as input for the lumpy script. I Manually entered in numbers to the lumpy script below for each sample but could easily save this output to a file and then parse file to obtain proper input to automate. 

# Extract the discordant paired-end alignments (discordants.bam)
samtools view -b -F 1294 sample.bam > sample.discordants.bam

# Extract the split-read alignments (splitters.bam)
samtools view -h sample.bam \
    | /usr/local/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.bam
    
#run lumpy
lumpy \
    -mw 4 \
    -tt 0 \
    -pe id:sample,bam_file:sample.discordants.bam,histo_file:sample.lib1.histo,mean:507,stdev:170,read_length:151,min_non_overlap:151,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:40 \
    -sr id:sample,bam_file:sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:40 \
    > /lumpy/sample_lumpy.vcf
```
gridss commands:
```
gridss.sh --jar /usr/local/gridss/scripts/gridss.jar --reference ref.fna --output sample_gridss.vcf.gz --labels sample --assembly sample_g.bam sample.bam
```
TODO: automated shell script to run all on desired samples 

### Merging Stuctural Variant Calls


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
- [ ] Obtain high quality variants for each strain 
- [ ] Extract all genes from reference using bioython 
- [ ] Write python code to modify reference gene sequences for each sample to contain ALT positions
- [ ] Concat all genes together to create consensus fasta sequence for each sample (header being sample name)
- [ ] Concat all samples together to create one multifastA 
- [ ] Use [MUSCLE](https://biopython.org/DIST/docs/api/Bio.Align.Applications._Muscle.MuscleCommandline-class.html) or [MAFFT](https://biopython.org/DIST/docs/api/Bio.Align.Applications._Mafft.MafftCommandline-class.html) python wrapper to create sequence alignment (phylip format) 
- [ ] Input sequence alignment into RAXML using following command line:

```
raxmlHPC -f a -x 43734 -p 89493 -# 100 -s snp_alignment.fasta -n phylo_result -m GTRGAMMA 
```
- [ ] Visualize tree using [FigTree](http://evomics.org/resources/software/molecular-evolution-software/figtree/) 

### Automating steps??? (get to this hopefully eventually)
- Snakemake pipeline one day 

