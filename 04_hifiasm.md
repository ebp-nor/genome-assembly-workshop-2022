# hifiasm tutorial

When choosing an assembler, you need to keep your data in mind. Since we want to create haplotype resolved assemblies, and we have both HiFi and Hi-C reads available, we are going to use **hifiasm**. Hifiasm is fast, easy to use, and creates high quality assemblies with longer contigs. To read more about how this software works, click [here.](https://github.com/chhylp123/hifiasm)

## Assembling with hifiasm

```
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --account=FIKS
#SBATCH --time=96:0:0
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=45000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

module --force purge

source /cluster/projects/path/to/conda.sh

eval "$(conda shell.bash hook)"

conda activate hifiasm

hifiasm -o $1 -t32  \
--h1 $2 \
--h2 $3 \
$4 \
1> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".err
```

```
awk '/^S/{print ">"$2"\n"$3}' in.gfa | fold > out.fa
```