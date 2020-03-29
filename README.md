Project Burk github Draft :-) 

#could use something like these to help brighten up readme and add color coding
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `Installation`
- ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) `Pipeline`
- ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) `Visualization`


### CF-BM-WGS-data-analysis


 - Repository containing pipeline and code involved in analyzing whole genome sequences of 19 Burkholderia multivorans isolates over 3 year time frame in cystic fibrosis patient. 

 ## Submodule descriptions 

- Submodules from master of each individual part of the analysis list contents/general part of project that submodule contains to keep everything organized.


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
- Manuscript
- Background info on CF/Burkholderia/bacterial evolution in vivo
- bacterial strain information/tables 
- methods for sequencing/isolation/culturing


