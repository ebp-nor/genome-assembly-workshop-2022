# GenomeScope2 tutorial

When creating de novo assemblies, there are a lot of considerations to take into account. What is the ploidy of the organism that you are assembling? What is the size of the genome? And what is the heterozygosity rate and repeat content? All these parameters, and more, can be determined by running **GenomeScope2**. For this tutorial we will be running the code below, but for more information about the software, you can click [*here*.](https://github.com/tbenavi1/genomescope2.0) 

## Creating a k-mer profile plot

```
#!/bin/bash
#SBATCH --job-name=genomescope
#SBATCH --account=ec146
#SBATCH --time=1:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate smudgescope

k=21
ploidy=2

mkdir -p tmp
echo $1 > FILES
[ -s reads.kmc_suf ] || kmc -k$k -t10 -m38 -ci1 -cs10000 @FILES reads tmp/

[ -s reads.histo ] || kmc_tools transform reads histogram reads.histo -cx10000

genomescope2 -i reads.histo -o output_ploidy1 -k $k -p 1 1> genomescope_ploidy1.out 2> genomescope_ploidy1.err
genomescope2 -i reads.histo -o output_ploidy2 -k $k -p 2 1> genomescope_ploidy2.out 2> genomescope_ploidy2.err
genomescope2 -i reads.histo -o output_ploidy4 -k $k -p 4 1> genomescope_ploidy4.out 2> genomescope_ploidy4.err

```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/genomescope`) with this content (with nano for instance):

sbatch /projects/ec146/scripts/run_genomescope.sh /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz
When you have done this, you can submit to the cluster by typing sh run.sh.

This should finish in a handful of minutes (when testing it ran for 2 minutes). You can monitor the progress with `squeue -u <username>.



## Interpreting your k-mer profile plot
