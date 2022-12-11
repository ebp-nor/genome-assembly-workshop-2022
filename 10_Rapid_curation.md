# Rapid curation tutorial

Although hifiasm is a great assembler, and YaHS can create chromosome length scaffolds, assembly errors do occur. Whether there are contigs that are misassembled, or scaffolds that are harder for the software to place, sometimes we need to manually curate the assemblies in order to reach our EBP-Nor assembly standards (you can read more about that here). To do this, we use the **Rapid curation** suite, developed by the GRIT-team at the Wellcome Sanger Institute, and **PretextView**, which you´ll learn more about in the last tutorial. If you want to read more about the code you´ll be using today, click [here](https://gitlab.com/wtsi-grit/rapid-curation/-/blob/main/README_software.md), and if you want to read more about why curation is so important for good quality reference genomes, click [here.](https://academic.oup.com/gigascience/article/10/1/giaa153/6072294) 

## Running the Rapid curation suite

### Run the suite

```
#!/bin/bash
#SBATCH --job-name=curation
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/fp/projects0101/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate curation

WORKDIR=$PWD/data
DESTDIR=$PWD/out
HICDIR=$2

export SINGULARITY_BIND="
$WORKDIR:/data,\
$HICDIR:/hic,\
$DESTDIR:/output,\
$TMP_DIR:/tmp
"

#hic
singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runHiC.sif -q 0 -s $1
rm $HOME/hic_done

#coverage
minimap2 -ax map-hifi \
         -t 5 data/ref.fa \
	$3 \
| samtools sort -@16 -O BAM -o coverage.bam

samtools view -b -F 256 coverage.bam > coverage_pri.bam

samtools index coverage_pri.bam

bamCoverage -b coverage_pri.bam -o coverage.bw

#gaps
singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runGap.sif -t $1

#repeats
singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runRepeat.sif -t $1  -s 10000

#telomers
#singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runTelo.sif -t $1 -s $3
#skipping telomers since they are not regular in budding yeast, at least not to our knowledge

#put it together
bigWigToBedGraph coverage.bw  /dev/stdout |PretextGraph -i out/out.pretext -n "PB coverage"

cat out/*_gap.bedgraph  | PretextGraph -i out/out.pretext -n "gaps"

#cat out/*_telomere.bedgraph |awk -v OFS="\t" '{$4 *= 1000; print}' | PretextGraph -i out/out.pretext -n "telomers"

bigWigToBedGraph  out/*_repeat_density.bw /dev/stdout | PretextGraph -i out/out.pretext -n "repeat density"
```


### Starting the script

```
mkdir -p data
mkdir -p out

cat  ../yahs/gsMetZobe_scaffolds_final.fa > data/ref.fa 

echo "/hic/hic_yeat.bam" > data/cram.fofn

sbatch /projects/ec146/scripts/run_rapidcuration.sh gsMetZobe /fp/projects01/ec146/data/genomic_data/hic/  /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz
```


### For information: converting fastq files to BAM
The rapid curation suite requires Hi-C reads to be in a BAM format. To create that, we did this:
```
#!/bin/bash
#SBATCH --job-name=convert_bam
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-node=10


module load picard/2.24.0-Java-11

java -Xms6g -Xmx6g -jar $EBROOTPICARD/picard.jar FastqToSam \
F1=ERR9503460_1_60x.fastq.gz \
F2=ERR9503460_1_60x.fastq.gz \
O=hic_yeast.bam \
SM=hic_yeast \
TMP_DIR=.
```

