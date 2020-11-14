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

# configfile: 'config.yaml'
# can be overwritten 
gtfs = {
    'rmsk': config.get('rmsk_gtf', '/home/schamori/data/snakepipes/GRCm38_98/annotation/rmsk.gtf'),
    'genes': config.get('genes_gtf', '/home/schamori/data/snakepipes/GRCm38_98/annotation/genes.gtf')
}

genes_bed = config.get('genes_bed', '/home/schamori/data/snakepipes/GRCm38_98/annotation/genes.bed')

rule all:
    input:
        'peakqc/plot.png',
        'peakqc/plot.svg',
        'heatmaps/genes.png',
        expand('heatmaps/{sample}_peaks.png', sample=SIGNAL_SAMPLES),
        # 'diffexp/table.csv',
        # 'diffexp/ma_plot.svg',
        expand('peaks/{sample}_{normalization}.{annot}.bed', normalization=['iggnormed', 'top1percent'], sample=SIGNAL_SAMPLES, annot=['genes', 'rmsk']),


include: 'rules/peakqc.smk'
include: 'rules/heatmap.smk'
include: 'rules/diffexp.smk'

rule scale_factor:
    input:
        ecoli="ecoli_mapping/filtered_bam/{sample}.filtered.bam"
    output:
        sf="scaleFactors/{sample}"
    conda:
        'env.yaml'
    shell: '''
        seqDepthDouble=$(samtools view -F 0x04 "{input.ecoli}" | wc -l)
        c=10000  # depends on the number of reads the ecoli mapping has.
        seqDepth=$((seqDepthDouble/2))
        if [[ "$seqDepth" -gt "1" ]]; then
            scale_factor=`echo "$c / $seqDepth" | bc -l`
            echo $scale_factor > {output.sf}
        else
            echo "$c" > {output.sf}
        fi
    '''

# generate scaled bedgraphs using ecoli-spikeins
rule bedgraph:
    input:
        mm10="mm10_mapping/filtered_bam/{sample}.filtered.bam",
        sf="scaleFactors/{sample}"
    output:
        bg="normalized/{sample}.bedgraph"
    conda:
        'env.yaml'
    log:
        "log/bedgraph/{sample}.log"
    shell: '''
        scale_factor=`cat {input.sf}`

        # no need to filter for unmapped reads (they are not counted here)
        bedtools genomecov -pc -bg -scale $scale_factor -ibam "{input.mm10}" > "{output.bg}" 2> {log}
    '''

# find peaks
rule seacr:
    input:
        ab="normalized/{AB}_{cond}_{repl}.bedgraph",
        igg="normalized/IgG_{cond}_{repl}.bedgraph"
    params:
        iggnormed=lambda wildcards, output: output['iggnormed'].replace('.bed', ''),
        top1percent=lambda wildcards, output: output['top1percent'].replace('.bed', '')
        # unpack(lambda wildcards, output: {key: value.replace('.bed', '') for key, value in output.items()})
    output:
        iggnormed="peaks/{AB}_{cond}_{repl}_iggnormed.bed",  # no need to exclude IgG in AB. It's excludeed before (AB_cond = ..)
        top1percent="peaks/{AB}_{cond}_{repl}_top1percent.bed"
    log:
        iggnormed="log/seacr/{AB}_{cond}_{repl}_iggnormed.log",
        top1percent="log/seacr/{AB}_{cond}_{repl}_top1percent.log"
    conda:
        'env.yaml'
    shell: '''
    SEACR_1.3.sh {input.ab} {input.igg} non stringent {params.iggnormed} {params.iggnormed} 2>&1 > {log.iggnormed}
        bedtools sort -i {params.iggnormed}.stringent.bed > {output.iggnormed}
        rm {params.iggnormed}.stringent.bed
    SEACR_1.3.sh {input.ab} 0.01 non stringent {params.top1percent} {params.top1percent} 2>&1 > {log.top1percent}
        bedtools sort -i {params.top1percent}.stringent.bed > {output.top1percent}
        rm {params.top1percent}.stringent.bed
    '''

rule generate_gene_only_gtf:
    'Generate GTF with only gene entries but containing TEs as well as protein_coding genes'
    input:
        lambda wildcards: gtfs[wildcards.annot]
    output:
        'ref/{annot}.gtf'
    # params:  # we could add custom filters for rmsk vs genes 
#         && $9 ~ /protein_coding/
    conda:
        'env.yaml'
    shell: '''
        awk -F $'\t' 'BEGIN {{ OFS=FS }} {{if ($3 == "gene" ) print $0;}}' {input} | sed 's/^chr//' | sort -k1,1 -k4,4n -s > {output}
    '''

rule bedtools_closest:
    input:
        peaks='peaks/{sample}.bed',  # both need to be sorted in the same manner.
        genes='ref/{annot}.gtf'
    output:
        'peaks/{sample}.{annot}.bed'
    params:
        field_offset="6"
    log:
        "log/bedtools_closest_{sample}.{annot}.log"
    conda:
        'env.yaml'
    shell: '''
         # adapted from snakepipes
         bedtools closest -d -a {input.peaks} -b {input.genes} \
         | cut -f1-{params.field_offset},$(( {params.field_offset} + 9 )),$(( {params.field_offset} + 10 )) \
         | awk -F $'\\t' '{{OFS=FS}} {{$7=gensub(".*gene_id \\"([^\\"]+)\\".*gene_name \\"([^\\"]+)\\".*", "\\\\1\\t\\\\2", $7) ; print}}' \
         > {output}  2> {log}
    '''
