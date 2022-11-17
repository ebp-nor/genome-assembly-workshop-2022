# Smudgeplot tutorial

The creators of GenomeScope2 created another way to visualize and estimate the ploidy and genome structure, using heterozygous k-mer pairs instead of k-mer frequency distribution. This software, called **Smudgeplot**, creates a “heatgraph” where a gradient from blue to yellow indicates the relative frequency of k-mer pairs. If the read coverage is good enough, you will have clear “smudges” of yellow which indicates the ploidy of your sequenced organism. The tutorial on how to create a smudgeplot can be found below, but if you want to read more about the software, you can click [here.](https://github.com/KamilSJaron/smudgeplot) 

## Creating a smudgeplot

```
#!/bin/bash
#SBATCH --job-name=smudgeplot
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=60G
#SBATCH --ntasks-per-node=10


#to get correct samtools: conda install -c bioconda samtools openssl=1.0

source /cluster/projects/path/to/conda.sh

eval "$(conda shell.bash hook)"

conda activate smudgeplot

#reads as fastq
#reads=$1
k=21
ploidy=2


mkdir -p tmp
ls *.fastq* > FILES
[ -s reads.kmc_suf ] || kmc -k$k -t10 -m38 -ci1 -cs10000 @FILES reads tmp/

[ -s reads.histo ] || kmc_tools transform reads histogram reads.histo -cx10000

L=$(smudgeplot.py cutoff reads.histo L)
U=$(smudgeplot.py cutoff reads.histo U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000


kmc_tools transform reads -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```

## Interpreting your smudgeplot