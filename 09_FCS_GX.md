# FCS-GX tutorial

Contaminants can end up in your assemblies in various different ways. Maybe someone touched samples without gloves? Maybe there were symbionts living on the organism when it was sampled? Or maybe the sample was contaminated during the sequencing run? Luckily, there are several genomic decontamination tools available, and the one we use in EBP-Nor is the **NCBI Foreign Contamination Screen (FCS)** tool suite (click [here](https://github.com/ncbi/fcs) to read more). This program suite can identify and remove contaminant sequences, whether it is adaptor sequences, vector contamination or foreign organisms. In today´s workshop, we are going to focus on the latter, and below you can find the code to run your own decontamination script. 

## Decontaminating the yeast assemblies

```
#!/bin/bash
#SBATCH --job-name=fcsgx
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-node=10

export SHM_LOC=/fp/projects01/ec146/opt/fcs

echo "GX_NUM_CORES=10" > env.txt

python3 /fp/projects01/ec146/opt/fcs/run_fcsgx.py --fasta $1 \
--gx-db  "${SHM_LOC}/gxdb/all" --split-fasta --tax-id $2 \
--gx-db-disk "${SHM_LOC}/gxdb/all.gxi" \
--container-engine singularity --image /fp/projects01/ec146/opt/fcs/fcsgx.sif
```

As we have done earlier, we have set up this script for you. Create a run.sh in your working folder (`/projects/ec146/work/<username>/fcsgx`) with this content (with `nano` for instance):

```
sbatch /projects/ec146/scripts/run_gcsgx.sh assembly.fasta \
taxonomy_id
```
You have to modify the run.sh script based on your assembly file and you have to find the taxonomy ID for *Metschnikowia zobellii* and input that.

Unfortunately this program requires a lot of memory to run (["approximately 470 GiB"](https://github.com/ncbi/fcs/wiki/FCS-GX)). If it is given unsufficient memory, the running time can increase by a factor of 10000x. On Fox, there are not that [many nodes](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/system-overview.md) with a lot of memory. The normal nodes have 501 GiB RAM, while the GPU accelerated nodes have up to 1006 GiB. Ideally, the job should have been allocated a bit more memory than what it strictly needs, but that is not easy here. Luckily, it should run in a handful of minutes if configured properly. 


After running the decontamination script, which foreign contaminants did you find?
