# Welcome to the Genome assembly, curation and validation workshop!

Insert introductary text.

## Our dataset

*Metschnikowia zobellii* is a yeast found in arctic coastal climates. It was first discovered as a water flea parasite, but have since been found in a wide range of organisms, including plants and other arthropods. *Metschnikowia zobellii* has a small genome, with only five chromosomes.  


**Why do we use a combination of HiFi and Hi-C reads?**

HiFi sequencing creates highly accurate circularized consensus reads. How are these reads generated? By ligating hairpin adapters, the DNA fragment that is being sequenced becomes a circle. This means that the machine can du multiple passes over the same DNA-sequence, to weed out any misread nucleotides. This is how HiFi reads can be so long, while remaining over 99.9% accurate. 

Hi-C sequencing is done to capture how the chromatin is folded within the cell nucleus. By ligating the folded DNA-strands, we can capture which loci are found in close proximity, and thus which parts of the DNA are found within the same chromosomes.

When combining these two, we can create haplotype resolved assemblies, meaning we can separate reads by maternal and paternal origin, without having access to parental data. In diploid, or polyploid organisms, this adds another level of information, and creates more accurate assemblies than a primary and alternate assembly would. 

## Package management

Write about conda, singularity and the modules we are using.

To load conda, do this:
```
eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 
```
## Infrastructure

We are using Educloud. Everyone need an account. Our project is ec146.

Write about the infrastructure that we are using in the workshop.
