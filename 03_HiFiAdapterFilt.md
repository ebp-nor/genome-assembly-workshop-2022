# HiFiAdapterFilt tutorial

Now that you have learned a bit more about your dataset you are *almost* ready to start the assembly process. Before you can start, you need to remove any remaining adapter sequences that may still be attached to your HiFi reads. To do this, we in EBP-Nor use **HiFiAdapterFilt**. To learn more about how this software works, click [here](https://github.com/sheinasim/HiFiAdapterFilt), and to do it yourself, follow the tutorial below.

## Filtering adapter sequences with HiFiAdapterFilt

```
#!/bin/bash
#SBATCH --job-name=hifiadaptfilt
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiadapterfilt

export PATH=/fp/projects01/ec146/opt/HiFiAdapterFilt/:$PATH
export PATH=/fp/projects01/ec146/opt/HiFiAdapterFilt/DB:$PATH

pbadapterfilt.sh -t 5
```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/hifiadaptfilt`) with this content (with nano for instance):
```
ln -s /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz . 
sbatch /projects/ec146/scripts/run_hifiadaptfilt.sh
```  
When you have done this, you can submit to the cluster by typing sh run.sh.

This should finish in a handful of minutes (when testing it ran for 1.5 minutes). You can monitor the progress with squeue -u <username>.

Smudgeplot produces several files in addition to the plot itself. You can for instance look at smudgeplot_verbose_summary.txt which contain the same information as the plot, but in text.



