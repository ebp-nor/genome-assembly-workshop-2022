# Merqury tutorial

Another way to validate your assemblies is by using **Merqury**. This k-mer based tool compares k-mers from the unassembled reads to the assemblies you have created, to find the degree of k-mer completeness, i.e. the percentage of the k-mers from the unassembled reads that are found within your assemblies. The advantage of this genome validation tool is that it can be used without any references, which is optimal when assessing de novo assemblies. If you want to learn more about Merqury, click [*here.*](https://github.com/marbl/merqury)

## Running Merqury

To run Merqury, create a new directory named `Merqury` and copy the code in the chunk below into a `run.sh` file:

```
#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
##SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=4500M
#SBATCH --ntasks-per-node=10


source /cluster/projects/nn8013k/programs/miniconda3/etc/profile.d/conda.sh

eval "$(conda shell.bash hook)"

conda activate merqury

#https://github.com/marbl/merqury

#$1 meryl db, $2 first asm $3 second asm, $4 output prefix 



mkdir -p $4
cd $4

merqury.sh $1 ../$2 ../$3 $4 > $4_merqury.out 2> $4_merqury.err
```

To run the script, write:

```
sbatch run.sh meryl.db first.fasta second.fasta prefix
```

## Interpreting a Merqury assembly spectrum plot