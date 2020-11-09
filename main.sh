#!/bin/bash
# inspired by https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye.html

CHROM_SIZES=mm10.chrom_sizes

# mouse mapping
DNA-mapping --local --fastqc --configFile ../config.yaml --input-dir fastq --output-dir mm10_mapping GRCm38_release98

# ecoli spikein mapping
DNA-mapping --local --fastqc --configFile ../config.yaml --input-dir fastq --output-dir ecoli_mapping --alignerOpts '--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10' ecoli_K12

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
snakemake --use-conda -j 30 -s $DIR/Snakefile
