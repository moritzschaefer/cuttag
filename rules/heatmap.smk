rule coordinate_sort_bam:
    input:
        "mm10_mapping/filtered_bam/{sample}.bam",
    output:
        "mm10_mapping/coordinate_sorted/{sample}.bam",
    conda:
        '../env.yaml'
    shell: '''
        samtools sort -o {output} {input}
    '''

rule index_bam:
    input:
        "mm10_mapping/coordinate_sorted/{sample}.bam",
    output:
        "mm10_mapping/coordinate_sorted/{sample}.bam.bai",
    conda:
        '../env.yaml'
    shell: '''
        samtools index {input}
    '''

rule generate_bw:
    input:
        bam="mm10_mapping/coordinate_sorted/{sample}.bam",
        index="mm10_mapping/coordinate_sorted/{sample}.bam.bai"
    output:
        "raw_coverage/{sample}.bw"
    conda:
        '../env.yaml'
    shell: '''
        bamCoverage -b {input.bam} -o {output}
    '''

rule generate_gene_matrix:
    input:
        bws=expand('raw_coverage/{sample}.bw', sample=SIGNAL_SAMPLES),
        gene_bed=config['heatmap_bed']
    output:
        'heatmaps/genes.mat.gz'
    conda:
        '../env.yaml'
    params:
        cores=20
    shell: '''
        computeMatrix scale-regions -S {input.bws} \
            -R {input.gene_bed} \
            --beforeRegionStartLength 3000 \
            --regionBodyLength 5000 \
            --afterRegionStartLength 3000 \
            --skipZeros -o {output} -p {params.cores}
    '''

rule plot_gene_heatmap:
    input:
        'heatmaps/genes.mat.gz'
    output:
        'heatmaps/genes.png'
    conda:
        '../env.yaml'
    shell: '''
        plotHeatmap -m {input} -out {output} --sortUsing sum
    '''

rule filter_peaks_summit:
    input:
        'peaks/{sample}_iggnormed.bed'
    output:
        'peaks/{sample}_summit.bed'
    conda:
        '../env.yaml'
    shell: '''
        awk '{{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}}' {input} \
            > {output}
    '''

rule generate_peak_matrix:
    input:
        raw_coverage="raw_coverage/{sample}.bw",
        peak_summit_bed="peaks/{sample}_summit.bed"
    output:
        'heatmaps/{sample}_peaks.mat.gz'
    params:
        cores=20
    conda:
        '../env.yaml'
    shell: '''
        computeMatrix reference-point -S "{input.raw_coverage}" \
                    -R {input.peak_summit_bed} \
                    --skipZeros -o {output} -p {params.cores} -a 3000 -b 3000 --referencePoint center
    '''

rule plot_peak_heatmap:
    input:
        'heatmaps/{sample}_peaks.mat.gz'
    output:
        'heatmaps/{sample}_peaks.png'
    conda:
        '../env.yaml'
    shell: '''
        plotHeatmap -m {input} -out {output} --sortUsing sum --startLabel "Peak Start" \
                    --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "{wildcards.sample}" 
    '''
# "${histName} ${repName}"
