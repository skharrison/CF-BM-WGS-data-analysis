
# Burkholderia WGS Analysis 

 - Repository containing all pieces and code involved in analyzing whole genome sequences of 19 Burkholderia multivorans isolates over 3 year time frame from sputum samples of a patient with cystic fibrosis.  

 ## Software Dependacies 

- [http://www.usadellab.org/cms/?page=trimmomatic][Trimmomatic] 


**Variant calling pipeline (software used to get vcf)**
- required software (QC, assembly, etc)
- python/perl code for allele frequency (AF) calculation, sequence error cutoff
- shell scripts used to execute all software and python (that way all command line options are displayed)

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


