How to use this pipeline:

- create a project folder
- Inside the project folder create a fastq folder
- Place your fastq's with the following naming schema
  {AB}_{condition}_{repl}_R[1|2].fastq
  - {AB} being the "IgG" or your AB-target of choice. Make sure it does NOT contain underscores
  - make sure condition does NOT contain underscores (i.e. use AGO2KO instead of AGO2_KO)
  - repl is the number of the biological replicate (e.g. 1, 2)

After placing the files, simply run

  conda activate snakePipes
  main.sh

from within the project directory (i.e. the directory that contains the fastq/ folder)

For the generation of the gene_heatmap, a bed file as generated by the snakepipes pipeline seems to work very well (in annotation/genes.bed).

Note: the ecoli-dna-mapping pipeline fails if there are 0 reads for a sample (e.g. IgGs usually have very few). However, it doesn't matter. Just make sure the scale_factors (they are stored in files) are to your liking.

Note: IgG is not supported anymore

# TODOs:
- adjust bowtie-params for multi-map sequencing (for TE analysis)
- integrate running snakepipes in the Snakefile
- make deseq analysis more flexible
- broad vs narrow peaks (can SEACR do narrow peaks)

