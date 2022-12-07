# YaHS tutorial

Congratulations, you have created your yeast assembly! But now you have to combine your contigs into scaffolds, and for that we use **YaHS**. YaHS stands for “yet another Hi-C scaffolding tool”, and as the name implies, there are a lot of Hi-C scaffolders out there. However, we in EBP-Nor choose to use YaHS because it is fast, creates more contiguous scaffolds, with better genome statistics compared to other widely used scaffolders. To learn more about this, click [here](https://github.com/c-zhou/yahs), otherwise scroll down to start scaffolding your haplotype resolved assemblies.

## Scaffolding with YaHS

```
#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate yahs

REF=$1

[ -s $REF.bwt ] || bwa index $REF

SAMPLE=$2

mkdir -p outs

[ -s hic_markdup.sort_n.bam ] || bwa mem -t 8 -R '@RG\tSM:$SAMPLE\tID:$SAMPLE' -5SPM $REF \
$3 $4 \
|samtools view -buS - |samtools sort -@1 -n -T tmp_n -O bam - \
|samtools fixmate -mr - -|samtools sort -@1 -T hic_tmp -O bam - |samtools markdup -rsS - -  2> hic_markdup.stats |samtools sort -n -@1 -n -T temp_n -O bam\
> hic_markdup.sort_n.bam

[ -s $REF.fai ] ||samtools faidx $REF

if [ -s $SAMPLE.bin ]; then
        yahs $REF $SAMPLE.bin -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
else
        yahs $REF hic_markdup.sort_n.bam -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
fi

```
