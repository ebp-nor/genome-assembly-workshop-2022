# GenomeScope2 tutorial

When creating de novo assemblies, there are a lot of considerations to take into account. What is the ploidy of the organism that you are assembling? What is the size of the genome? And what is the heterozygosity rate and repeat content? All these parameters, and more, can be determined by running **GenomeScope2**. For this tutorial we will be running the code below, but for more information about the software, you can click [*here*.](https://github.com/tbenavi1/genomescope2.0) 

## Creating a k-mer profile plot

```
#!/bin/bash
#SBATCH --job-name=genomescope
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=4500M
#SBATCH --ntasks-per-node=10


#to get correct samtools: conda install -c bioconda samtools openssl=1.0

source /cluster/projects/path/to/conda.sh

eval "$(conda shell.bash hook)"

conda activate genomescope

k=21
ploidy=2

mkdir -p tmp
ls *.fastq* > FILES
[ -s reads.kmc_suf ] || kmc -k$k -t10 -m38 -ci1 -cs10000 @FILES reads tmp/

[ -s reads.histo ] || kmc_tools transform reads histogram reads.histo -cx10000

genomescope2 -i reads.histo -o output_ploidy1 -k $k -p 1 1> genomescope_ploidy1.out 2> genomescope_ploidy1.err
genomescope2 -i reads.histo -o output_ploidy2 -k $k -p 2 1> genomescope_ploidy2.out 2> genomescope_ploidy2.err
genomescope2 -i reads.histo -o output_ploidy4 -k $k -p 4 1> genomescope_ploidy4.out 2> genomescope_ploidy4.err

```


## Interpreting your k-mer profile plot