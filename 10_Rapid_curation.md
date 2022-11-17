# Rapid curation tutorial

Insert Rapid curation tutorial from Notion.

## Running the Rapid curation suite

### Run HiC

```
#!/bin/bash
#SBATCH --job-name=run_hic
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=8000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10


WORKDIR=$PWD/data
DESTDIR=$PWD/out
HICDIR=$3                                                            
TMP_DIR=$2

export SINGULARITY_BIND="
$WORKDIR:/data,\
$HICDIR:/hic,\
$DESTDIR:/output,\
$TMP_DIR:/tmp
"

singularity run /cluster/projects/nn8013k/opt/rapid-curation/rapid_hic_software/runHiC.sif -q 0 -s $1
rm $HOME/hic_done
```

### Run coverage

```
#!/bin/bash
#SBATCH --job-name=minimap_mall
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=3900M
#SBATCH --ntasks-per-node=16



source /cluster/projects/nn8013k/programs/miniconda3/etc/profile.d/conda.sh

eval "$(conda shell.bash hook)"

conda activate minimap


minimap2 -ax map-hifi \
         -t 16 data/ref.fa \
	$1 \
| samtools sort -@16 -O BAM -o coverage.bam

samtools view -b -F 256 coverage.bam > coverage_pri.bam

samtools index coverage_pri.bam

bamCoverage -b coverage_pri.bam -o coverage.bw
```


### Run gap

```
#!/bin/bash
#SBATCH --job-name=run_gap
#SBATCH --account=FIKS
#SBATCH --time=1:0:0
#SBATCH --mem-per-cpu=8000M
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-node=10


WORKDIR=$PWD/data
DESTDIR=$PWD/out
TMP_DIR=$2


export SINGULARITY_BIND="
$WORKDIR:/data,\
$DESTDIR:/output,\
$TMP_DIR:/tmp,\
"

singularity run /cluster/projects/nn8013k/opt/rapid-curation/rapid_hic_software/runGap.sif -t $1
```


### Run repeat

```
#!/bin/bash
#SBATCH --job-name=run_repeat
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
##SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=8000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

WORKDIR=$PWD/data
DESTDIR=$PWD/out
TMP_DIR=$2  

export SINGULARITY_BIND="
$WORKDIR:/data,\
$DESTDIR:/output,\
$TMP_DIR:/tmp,\
"

singularity run /cluster/projects/nn8013k/opt/rapid-curation/rapid_hic_software/runRepeat.sif -t $1  -s 10000
```


### Run telomere

```
#!/bin/bash
#SBATCH --job-name=run_telo
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=8000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

WORKDIR=$PWD/data
DESTDIR=$PWD/out
TMP_DIR=$2


export SINGULARITY_BIND="
$WORKDIR:/data,\
$DESTDIR:/output,\
$TMP_DIR:/tmp,\
"

singularity run /cluster/projects/nn8013k/opt/rapid-curation/rapid_hic_software/runTelo.sif -t $1 -s $3
```


### Running all the scripts at once

```
mkdir -p data
mkdir -p out

cat kcLamFluv2.h1.decon.fasta > data/ref.fa 

ls /cluster/projects/nn8013k/results/species/Lampetra_fluviatilis/kcLamFluv2/genomic_data/hic/Sample_Omni-C-RiverLamprey/*bam |sed "s|/cluster/projects/nn8013k/results/species/Lampetra_fluviatilis/kcLamFluv2/genomic_data/hic/Sample_Omni-C-RiverLamprey|/hic/|g" > data/cram.fofn


#sbatch /cluster/projects/nn8013k/scripts/curation/run_runhic.sh kcLamFluv2_h1 /cluster/work/users/benedga/tmp /cluster/projects/nn8013k/results/species/Lampetra_fluviatilis/kcLamFluv2/genomic_data/hic/Sample_Omni-C-RiverLamprey
#sbatch /cluster/projects/nn8013k/scripts/curation/run_coverage.sh  /cluster/projects/nn8013k/results/species/Lampetra_fluviatilis/kcLamFluv2/genomic_data/pacbio/hifiadapterfilt/concat.filt.fastq.gz
#sbatch /cluster/projects/nn8013k/scripts/curation/run_rungap.sh kcLamFluv2_h1 /cluster/work/users/benedga/tmp
#sbatch /cluster/projects/nn8013k/scripts/curation/run_runrepeat.sh kcLamFluv2_h1 /cluster/work/users/benedga/tmp
#sbatch /cluster/projects/nn8013k/scripts/curation/run_runtelo.sh kcLamFluv2_h1 /cluster/work/users/benedga/tmp TTAGGG
```


### Preparing for PretextView

```
source  /cluster/projects/nn8013k/programs/miniconda3/etc/profile.d/conda.sh

eval "$(conda shell.bash hook)"

conda activate pretext

bigWigToBedGraph coverage.bw  /dev/stdout |PretextGraph -i out/out.pretext -n "PB coverage"

cat out/*_gap.bedgraph  | PretextGraph -i out/out.pretext -n "gaps"

cat out/*_telomere.bedgraph |awk -v OFS="\t" '{$4 *= 1000; print}' | PretextGraph -i out/out.pretext -n "telomers"

bigWigToBedGraph  out/*_repeat_density.bw /dev/stdout | PretextGraph -i out/out.pretext -n "repeat density"
```