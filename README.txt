How to use this pipeline:

- create a project folder
- Inside the project folder create a fastq folder
- Place your fastq's with the following naming schema
  {AB}_{condition}_{repl}_R[1|2].fastq
  - {AB} being the "IgG" or your AB-target of choice. Make sure it does NOT contain underscores
  - make sure condition does NOT contain underscores (i.e. use AGO2KO instead of AGO2_KO)
  - repl is the number of the biological replicate (e.g. 1, 2)

After placing the files, simply run
conda activate snakePipes3
snakemake --use-conda -j 25 --snakefile ../Snakefile
from within the project directory

# TODOs:
- adjust bowtie-params for multi-map sequencing (for TE analysis)
- integrate running snakepipes in the Snakefile
- compare their data run/pipeline!
- heatmaps (adjust gene_bed in config since it contains >100.000 transcripts at the moment...)
- make deseq analysis more flexible

