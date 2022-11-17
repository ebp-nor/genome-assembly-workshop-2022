# FCS_GX tutorial

Contaminants can end up in your assemblies in various different ways. Maybe someone touched samples without gloves? Maybe there were symbionts living on the organism when it was sampled? Or maybe the sample was contaminated during the sequencing run? Luckily, there are several genomic decontamination tools available, and the one we use in EBP-Nor is the **NCBI Foreign Contamination Screen (FCS)** tool suite (click [here](https://github.com/ncbi/fcs) to read more). This program suite can identify and remove contaminant sequences, whether it is adaptor sequences, vector contamination or foreign organisms. In today´s workshop, we are going to focus on the latter, and below you can find the code to run your own decontamination script. 

## Decontaminating the yeast assemblies

After running the decontamination script, which foreign contaminants did you find?