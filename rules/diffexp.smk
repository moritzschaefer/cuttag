import re

rule deseq:
    input:
        peaks=expand('peaks/{sample}.stringent.bed', sample=SIGNAL_SAMPLES),
        bams=expand('mm10_mapping/filtered_bam/{sample}.bam', sample=SIGNAL_SAMPLES)
    output:
        ''
    params:
        sample_names=SIGNAL_SAMPLES,
        condition_names=[re.sub('_\d+$', '', s) for s in SIGNAL_SAMPLES]
    conda:
        '../env.yaml'
    script:
        '../scripts/deseq.R'
