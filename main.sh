#!/bin/bash
# inspired by https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye.html

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

CHROM_SIZES=mm10.chrom_sizes

# mouse mapping
# aligner-opts defined in config.yaml
DNA-mapping --local --configFile $DIR/config.yaml --input-dir fastq --output-dir mm10_mapping GRCm38_98

# ecoli spikein mapping
# the protocol is ambiguous for this analysis
DNA-mapping --local --configFile $DIR/config.yaml --input-dir fastq --output-dir ecoli_mapping --alignerOpts '--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10' ecoli_K12

snakemake --use-conda -j 30 -s $DIR/Snakefile
