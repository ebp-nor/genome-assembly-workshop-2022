# Welcome to the Genome assembly, curation and validation workshop!

Today you´ll learn how to assemble a whole genome, EBP-Nor style! 

After attending the workshop you should:
- know about the most-used approaches for genome assembly
- be able to assess information inherit in sequencing reads
- be able to validate genome assemblies
- know about manual curation of assemblies


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

For the different analyses we are doing, we will use [Educloud](https://www.uio.no/english/services/it/research/platforms/edu-research/). To use it, you need an account which you can get here: [https://research.educloud.no/register](https://research.educloud.no/register). The project we are using in this course is ec146, so please ask for access to that one, and we will let you in. 

We will do the work in this course at `/projects/ec146/work` on [Fox](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/) which is the HPC part of Educloud. After creating an account, you can log in using `ssh <educloud-username>@fox.educloud.no`. You will be prompted for a One-Time Code for a 2-factor authenticator app (Microsoft Authenticator) and your Fox/Educloud password.

On Fox we will submit jobs/analyses as job scripts. This is for a system called SLURM. Basically, this are instructions to the system for what kind of analysis we are running, or more concretely, how much memory and computing power we need. 

A generic job script might look like this (copied from [https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/jobs/job-scripts.md](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/jobs/job-scripts.md):
```
#!/bin/bash

# Job name:
#SBATCH --job-name=YourJobname
#
# Project:
#SBATCH --account=ecXXX
#
# Wall time limit:
#SBATCH --time=DD-HH:MM:SS
#
# Other parameters:
#SBATCH ...

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load SomeProgram/SomeVersion
module list

## Do some work:
YourCommands
```


