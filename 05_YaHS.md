# YaHS tutorial

Congratulations, you have created your yeast assembly! But now you have to combine your contigs into scaffolds, and for that we use **YaHS**. YaHS stands for “yet another Hi-C scaffolding tool”, and as the name implies, there are a lot of Hi-C scaffolders out there. However, we in EBP-Nor choose to use YaHS because it is fast, creates more contiguous scaffolds, with better genome statistics compared to other widely used scaffolders. To learn more about this, click [here](https://github.com/c-zhou/yahs), otherwise scroll down to start scaffolding your haplotype resolved assemblies.

## Scaffolding with YaHS

```
#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --account=FIKS
#SBATCH --time=72:0:0
#SBATCH --mem-per-cpu=4500M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16


source /cluster/projects/path/to/conda.sh

eval "$(conda shell.bash hook)"

conda activate yahs

export TMPDIR=/cluster/path/to/tmp

PATH=/cluster/projects/nn8013k/opt/yahs:$PATH

REF=$1

[ -s $REF.bwt ] || bwa index $REF

SAMPLE=$2

mkdir -p outs

[ -s hic_markdup.sort_n.bam ] || bwa mem -t 10 -R '@RG\tSM:$SAMPLE\tID:$SAMPLE' -5SPM $REF \
$3 $4 \
|samtools view -buS - |samtools sort -@3 -n -T tmp_n -O bam - \
|samtools fixmate -mr - -|samtools sort -@3 -T hic_tmp -O bam - |samtools markdup -rsS - -  2> hic_markdup.stats |samtools sort -n -@3 -n -T temp_n -O bam\
> hic_markdup.sort_n.bam

[ -s $REF.fai ] ||samtools faidx $REF

if [ -s $SAMPLE.bin ]; then
        yahs $REF $SAMPLE.bin -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
else
        yahs $REF hic_markdup.sort_n.bam -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
fi


#yahs $REF hic_markdup.sort_n.bam  -o $SAMPLE \
#1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err

gfastats $SAMPLE_scaffolds_final.fa > $SAMPLE_scaffolds_final.stats

[ -s ${SAMPLE}_scaffolds_final.fa.fai ] ||samtools faidx ${SAMPLE}_scaffolds_final.fa

cut -f1-2 ${SAMPLE}_scaffolds_final.fa.fai >${SAMPLE}_scaffolds_final.chrom.sizes

(juicer_pre ${SAMPLE}.bin ${SAMPLE}_scaffolds_final.agp ${REF}.fai 2> tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T . --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' ${SAMPLE}_scaffolds_final.chrom.sizes; awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' alignments_sorted.txt) | PretextMap -o ${SAMPLE}.pretext
```