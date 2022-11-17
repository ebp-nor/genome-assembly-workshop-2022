# HiFiAdapterFilt tutorial

Now that you have learned a bit more about your dataset you are *almost* ready to start the assembly process. Before you can start, you need to remove any remaining adapter sequences that may still be attached to your HiFi reads. To do this, we in EBP-Nor use **HiFiAdapterFilt**. To learn more about how this software works, click [here](https://github.com/sheinasim/HiFiAdapterFilt), and to do it yourself, follow the tutorial below.

## Filtering adapter sequences with HiFiAdapterFilt

```
#!/bin/bash
#SBATCH --job-name=hifiadapterfilt
#SBATCH --account=FIKS
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=4500M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16


module load GCC/11.2.0 CMake/3.21.1-GCCcore-11.2.0

source /cluster/projects/path/to/conda.sh

eval "$(conda shell.bash hook)"
conda activate hifiadapterfilt

export PATH=/cluster/projects/nn9244k/olekto/projects/EBP-Nor/opt/HiFiAdapterFilt/:$PATH
export PATH=/cluster/projects/nn9244k/olekto/projects/EBP-Nor/opt/HiFiAdapterFilt/DB:$PATH
export TMPDIR=/cluster/work/users/olekto/tmp

pbadapterfilt.sh -t 16
```