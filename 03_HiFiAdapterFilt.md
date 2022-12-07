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

export PATH=/projects/ec146/opt/HiFiAdapterFilt/:$PATH
export PATH=/projects/ec146/opt/HiFiAdapterFilt/DB:$PATH

pbadapterfilt.sh -t 5
```
