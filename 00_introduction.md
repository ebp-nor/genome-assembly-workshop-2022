# Welcome to the Genome assembly, curation and validation workshop!

Insert introductary text.

## Our dataset

Write about the yeast


**Why do we use a combination of HiFi and Hi-C reads?**

HiFi sequencing creates circularized consensus reads. What this means is that by ligating hairpin adapters, the DNA fragment that is being sequenced by the Sequel II becomes a circle. This means that the machine can du multiple passes over the same DNA-sequence, to weed out any misread nucleotides. This is how HiFi reads can be so long, while remaining over 99.9% accurate. 

Hi-C sequencing is done to capture how the chromatin is folded within the cell nucleus. By ligating the folded DNA-strands, we can capture which loci are found in close proximity, and thus which parts of the DNA are found within the same chromosomes.

When combining these two, we can create haplotype resolved assemblies, meaning we can separate reads by maternal and paternal origin, without having access to parental data. In diploid, or polyploid organisms, this adds another level of information, and creates more accurate assemblies than a primary and alternate assembly would. 

## Package management

Write about conda, singularity and the modules we are using.

## Infrastructure

Write about the infrastructure that we are using in the workshop.