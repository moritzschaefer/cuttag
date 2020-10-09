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
AB_cond = list(set([re.match('[^_]+_[^_]+', f).group() for f in os.listdir('fastq') if 'IgG' not in f and 'fastq.gz' in f]))

rule all:
    input:
        'peakqc/plot.png',
        'peakqc/plot.svg',

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
            touch "normalized/$sample.bedgraph"  # generate empty file anyways
            echo "seqDepth=0 for $sample"
        fi
    '''

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

rule sort_bam:
    input:
        sample_bam="mm10_mapping/Bowtie2/{AB}_{cond}_{repl}.bam",
    output:
        sorted_bam="mm10_mapping/sorted_bam/{AB}_{cond}_{repl}.bam",
    conda:
        'env.yaml'
    shell: '''
        samtools sort -n "{input.sample_bam}" > "{output.sorted_bam}"
    '''
# get quality control values
rule peak_qc:
    input:
        sample="peaks/{AB}_{cond}_{repl}_{type}.stringent.bed",
        # other_replicate="peaks/{AB}_{cond}_{repl}_{type}.stringent.bed",
        other_replicate=lambda wildcards: f"peaks/{wildcards.AB}_{wildcards.cond}_{(int(wildcards.repl) % 2) + 1}_{wildcards.type}.stringent.bed",
        sample_bam_sorted="mm10_mapping/sorted_bam/{AB}_{cond}_{repl}.bam",
    output:
        "peakqc/{AB}_{cond}_{repl,\d+}_{type}.csv"
    conda:
        'env.yaml'
    shell: '''
        # replicate overlap
        overlap_peaks=$(bedtools intersect -wa -a "{input.sample}" -b "{input.other_replicate}" | wc -l)
        total_peaks=$(cat "{input.sample}" | wc -l)
        repl_percent=$(echo "100 * $overlap_peaks / $total_peaks" | bc -l)

        # TODO peak width

        # FRagment proportion in Peaks regions (FRiPs)
        mappedFragsDouble=$(samtools view -F 0x04 "{input.sample_bam_sorted}" | wc -l)
        mappedFrags=$((mappedFragsDouble/2))
        # weirdly, 'bedtools intersect' reports exactly twice as many reads as pairtobed. seems like bedtools intersect considers the whole fragment and always reports both reads per fragment
        mappedFragsInPeaks=$(bedtools pairtobed -type ospan -abam {input.sample_bam_sorted} -b {input.sample} -bedpe | wc -l)
        frip=$(echo "100 * $mappedFragsInPeaks / $mappedFrags" | bc -l)
        echo "antibody,condition,replicate,type,mapped_fragments,mapped_fragments_in_peaks,frips,total_peaks,overlap_peaks,repl_percent" > {output}
        echo "{wildcards.AB},{wildcards.cond},{wildcards.repl},{wildcards.type},$mappedFrags,$mappedFragsInPeaks,$frip,$total_peaks,$overlap_peaks,$repl_percent" >> {output}
    '''
rule plot_peak_qcs:
    input:
        csvs=expand('peakqc/{AB_cond}_{repl}_{type}.csv', AB_cond=AB_cond, repl=[1, 2], type=['iggnormed', 'top1percent']),
        beds=expand('peaks/{AB_cond}_{repl}_{type}.stringent.bed', AB_cond=AB_cond, repl=[1, 2], type=['iggnormed', 'top1percent'])
    output:
        png='peakqc/plot.png',
        svg='peakqc/plot.svg',
    conda:
        'env.yaml'
    script:
        'plot_peak_qc.py'
