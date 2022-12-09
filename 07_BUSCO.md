# BUSCO tutorial

How do you measure how “complete” your assembly is? Since we know our yeast assemblies are about the expected size (which we estimated using GenomeScope2), and that our assembly statistics are good (which we found out by using gfastats), we now want to know that the genes we expect to find are there. **BUSCO** is a tool that can be used to compare your assemblies to lists of near universal orthologous genes (i.e. genes with common ancestry and function), to find out if you have managed to assemble them correctly. The assumption is that if you use the correct lineage dataset (based on what organism you are assembling), you´ll be able to find most of the genes within your assembly. If you want to read more about BUSCO, click [*here.*](https://busco.ezlab.org/busco_userguide.html)


## Running BUSCO

This is the script we use to run BUSCO:

```
#!/bin/sh
#SBATCH --job-name=busco5
#SBATCH --partition=normal
#SBATCH --account=FIKS
#SBATCH --time=5:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G #Memory per node



module --force purge

source FIKS

eval "$(conda shell.bash hook)"

conda activate busco5


prefix=${1%.*}

ln -s $1 ${prefix}.fna

mkdir -p busco5_${2}_${prefix}
origdir=$PWD

cd busco5_${2}_${prefix}

echo $PWD

echo ${1}

busco -c 10 -i ${origdir}/$prefix.fna -l  /cluster/projects/FIKS/opt/busco_dbs/lineages/${2}_odb10 -o assembly -m genome  --offline > busco.out 2> busco.err

``` 

Create a new directory in your work area named `BUSCO`. Make a new `run.sh` file with `nano run.sh`, and copy the code below into that file:

```
ln -s path/to/assembly_hap1.fasta .
ln -s path/to/assembly_hap2.fasta .
sbatch /projects/ec146/scripts/run_busco.sh assembly_hap1.fasta
sbatch /projects/ec146/scripts/run_busco.sh assembly_hap2.fasta
```

When you have done this, you can submit to the cluster by typing `sh run.sh`.

Til Ole, jeg synes vi burde hardkode lineage inn i scriptet så de blir tvingt til å faktisk se på opsjonene, og skjønne hva -l spesifiserer.

## Reviewing the BUSCO results

Before you start reviewing your results, look at the [list of lineages](https://busco-data.ezlab.org/v5/data/lineages/) found within the OrthoDB database. Which lineage would you pick for the yeast assemblies? Is it the same as the one we specified in the script above?

**Try to answer these questions by reviewing your BUSCO results:**

1. Which of the haplotypes has the highest percentage of complete BUSCO genes?

2. What do you think it means for the assembly quality if the number of fragmented, missing or duplicated BUSCO genes are high? Discuss with the people on your table. 

