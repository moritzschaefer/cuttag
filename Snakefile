import os
import re

# could use this for running snakepipes first.. https://snakepipes.readthedocs.io/en/latest/content/advanced_usage.html#calling-snakemake-directly-using-the-snakefiles
# subworkflow otherworkflow:
    # workdir:
    #     "../path/to/otherworkflow"
    # snakefile:
    #     "../path/to/otherworkflow/Snakefile"
    # configfile:
    #     "path/to/custom_configfile.yaml"

# get antibody_condition combinations
ALL_SAMPLES = list(set([re.match('[^_]+_[^_]+_[^_]+', f).group() for f in os.listdir('fastq') if 'fastq.gz' in f]))
SIGNAL_SAMPLES = list(set([s for s in ALL_SAMPLES if 'IgG' not in s]))
AB_cond = list(set([re.match('[^_]+_[^_]+', s).group() for s in SIGNAL_SAMPLES]))

configfile: 'config.yaml'


rule all:
    input:
        'peakqc/plot.png',
        'peakqc/plot.svg',
        'heatmaps/genes.png',
        expand('heatmaps/{sample}_peaks.png', sample=SIGNAL_SAMPLES),


include: 'rules/peakqc.smk'
include: 'rules/heatmap.smk'
include: 'rules/diffexp.smk'

# generate scaled bedgraphs using ecoli-spikeins
rule bedgraph:
    input:
        mm10="mm10_mapping/Bowtie2/{sample}.bam",
        ecoli="ecoli_mapping/Bowtie2/{sample}.bam"
    output:
        bg="normalized/{sample}.bedgraph",
        sf="scaleFactors/{sample}"
    conda:
        'env.yaml'
    log:
        "log/bedgraph/{sample}.log"
    shell: '''
        seqDepthDouble=$(samtools view -F 0x04 "{input.ecoli}" | wc -l)
        seqDepth=$((seqDepthDouble/2))
        if [[ "$seqDepth" -gt "1" ]]; then  # only continue if we have spike_in reads
            scale_factor=`echo "10000 / $seqDepth" | bc -l`
            echo $scale_factor > {output.sf}

            # no need to filter for unmapped reads (they are not counted here)
            bedtools genomecov -pc -bg -scale $scale_factor -ibam "{input.mm10}" > "{output.bg}" 2> {log}
        else
            touch {output.bg}  # generate empty file anyways
            echo "-1" > {output.sf}
            echo "seqDepth=0 for {wildcards.sample}" > {log}
        fi
    '''

# find peaks
rule seacr:
    input:
        ab="normalized/{AB}_{cond}_{repl}.bedgraph",
        igg="normalized/IgG_{cond}_{repl}.bedgraph"
    params:
        iggnormed=lambda wildcards, output: output['iggnormed'].replace('.stringent.bed', ''),
        top1percent=lambda wildcards, output: output['top1percent'].replace('.stringent.bed', '')
        # unpack(lambda wildcards, output: {key: value.replace('.stringent.bed', '') for key, value in output.items()})
    output:
        iggnormed="peaks/{AB}_{cond}_{repl}_iggnormed.stringent.bed",  # no need to exclude IgG in AB. It's excludeed before (AB_cond = ..)
        top1percent="peaks/{AB}_{cond}_{repl}_top1percent.stringent.bed"
    log:
        iggnormed="log/seacr/{AB}_{cond}_{repl}_iggnormed.log",
        top1percent="log/seacr/{AB}_{cond}_{repl}_top1percent.log"
    conda:
        'env.yaml'
    shell: '''
        SEACR_1.3.sh {input.ab} {input.igg} non stringent {params.iggnormed} 2>&1 > {log.iggnormed}
        SEACR_1.3.sh {input.ab} 0.01 non stringent {params.top1percent} 2>&1 > {log.top1percent}
        '''
