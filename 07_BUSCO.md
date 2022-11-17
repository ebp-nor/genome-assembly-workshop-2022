# BUSCO tutorial

How do you measure how “complete” your assembly is? Since we know our yeast assemblies are about the expected size (which we estimated using GenomeScope2), and that our assembly statistics are good (which we found out by using gfastats), we now want to know that the genes we expect to find are there. **BUSCO** is a tool that can be used to compare your assemblies to lists of near universal orthologous genes (i.e. genes with common ancestry and function), to find out if you have managed to assemble them correctly. The assumption is that if you use the correct lineage dataset (based on what organism you are assembling), you´ll be able to find most of the genes within your assembly. If you want to read more about BUSCO, click [here.](https://busco.ezlab.org/busco_userguide.html)

## Reviewing the BUSCO results

Before you start reviewing your results, look at the [list of lineages](https://busco-data.ezlab.org/v5/data/lineages/) found within the OrthoDB database. Which lineage would you pick for the yeast assemblies? 