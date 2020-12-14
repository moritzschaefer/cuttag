import re


rule deseq_init:
    '''
    Combine peaks and run DESeq analysis
    '''
    input:
        peaks=expand('peaks/{sample}_top{{percent}}percent.bed', sample=SIGNAL_SAMPLES),  # TODO for now we use top1percent...
        bams=expand('mm10_mapping/filtered_bam/{sample}.filtered.bam', sample=SIGNAL_SAMPLES)
    output:
        rds="diffexp/all_top{percent}percent.rds",
        master_peaks="diffexp/master_peaks_top{percent}percent.bed"
    params:
        sample_names=SIGNAL_SAMPLES,
        condition_names=CONDITION_NAMES,
        batch_numbers=[re.search('_(\d+)$', s).groups()[0] for s in SIGNAL_SAMPLES]
    conda:
        '../env.yaml'
    script:
        '../scripts/deseq_init.R'

rule deseq_results:
    '''
    Get and store results. Output table contains the following values (from the henipipe protocol):
    "
    - First 4 columns: raw reads counts after filtering the peak regions with low counts
    - Second 4 columns: normalized read counts eliminating library size difference.
    - Remaining columns: differential detection results.
    "

    '''
    input:
        rds='diffexp/all_top{percent}percent.rds'
    params:
        control=NORMALIZER
    output:
        ma_plot='diffexp/ma_plot_{target}_top{percent}percent.svg',
        table='diffexp/table_{target}_top{percent}percent.csv'
    conda:
        '../env.yaml'
    script:
        '../scripts/deseq_results.R'
