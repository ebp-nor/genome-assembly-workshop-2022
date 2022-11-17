# FCS_GX tutorial

Contaminants can end up in your assemblies in various different ways. Maybe someone touched samples without gloves? Maybe there were symbionts living on the organism when it was sampled? Or maybe the sample was contaminated during the sequencing run? Luckily, there are several genomic decontamination tools available, and the one we use in EBP-Nor is the **NCBI Foreign Contamination Screen (FCS)** tool suite (click [here](https://github.com/ncbi/fcs) to read more). This program suite can identify and remove contaminant sequences, whether it is adaptor sequences, vector contamination or foreign organisms. In today´s workshop, we are going to focus on the latter, and below you can find the code to run your own decontamination script. 

## Decontaminating the yeast assemblies

```
#!/bin/bash
#SBATCH --job-name=fcsgx
#SBATCH --account=nn8013k
#SBATCH --time=5:0:0
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=60G 
#SBATCH --cpus-per-task=10

export SHM_LOC=/cluster/projects/nn8013k/opt/fcs

echo "GX_NUM_CORES=10" > env.txt

python3 /cluster/projects/nn8013k/opt/fcs/dist/run_fcsgx.py --fasta $1 --out-dir ./gx_out/ \
--gx-db  "${SHM_LOC}/gxdb/all" --split-fasta --tax-id $2 \
--gx-db-disk "${SHM_LOC}/gxdb/all.gxi" \
--container-engine=singularity --image=/cluster/projects/nn8013k/opt/fcs/fcsgx.sif
```


After running the decontamination script, which foreign contaminants did you find?