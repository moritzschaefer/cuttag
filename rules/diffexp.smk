import re

rule deseq:
    input:
        peaks=expand('peaks/{sample}_iggnormed.stringent.bed', sample=SIGNAL_SAMPLES),
        bams=expand('mm10_mapping/filtered_bam/{sample}.filtered.bam', sample=SIGNAL_SAMPLES)
    output:
        ma_plot='diffexp/ma_plot.svg',
        table='diffexp/table.csv'
    params:
        sample_names=SIGNAL_SAMPLES,
        condition_names=[re.sub('_\d+$', '', s) for s in SIGNAL_SAMPLES]
    conda:
        '../env.yaml'
    script:
        '../scripts/deseq.R'
