# hifiasm tutorial

When choosing an assembler, you need to keep your data in mind. Since we want to create haplotype resolved assemblies, and we have both HiFi and Hi-C reads available, we are going to use **hifiasm**. Hifiasm is fast, easy to use, and creates high quality assemblies with longer contigs. To read more about how this software works, click [here.](https://github.com/chhylp123/hifiasm)

## Assembling with hifiasm

```
#!/bin/bash
#SBATCH --job-name=smudgeplot
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiasm

hifiasm -o $1 -t5  \
--h1 $2 \
--h2 $3 \
$4 \
1> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".err
```

```
awk '/^S/{print ">"$2"\n"$3}' in.gfa | fold > out.fa
```


We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/smudgeplot`) with this content (with nano for instance): 
 
```
sbatch /projects/ec146/scripts/run_smudgeplot.sh /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz
```

When you have done this, you can submit to the cluster by typing `sh run.sh`.
 
This should finish in a handful of minutes (when testing it ran for 1.5 minutes). You can monitor the progress with `squeue -u <username>`.
