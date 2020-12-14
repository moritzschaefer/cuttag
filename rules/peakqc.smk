rule name_sort_bam:
    input:
        sample_bam="mm10_mapping/filtered_bam/{AB}_{cond}_{repl}.filtered.bam",
    output:
        sorted_bam="mm10_mapping/name_sorted/{AB}_{cond}_{repl}.bam",
    conda:
        '../env.yaml'
    shell: '''
    samtools sort -n "{input.sample_bam}" > "{output.sorted_bam}"
    '''


# get quality control values
rule peak_qc:
    input:
        sample="peaks/{AB}_{cond}_{repl}_top{percent}percent.bed",
        # other_replicate="peaks/{AB}_{cond}_{repl}_{type}.bed",
        other_replicate=lambda wildcards: f"peaks/{wildcards.AB}_{wildcards.cond}_{(int(wildcards.repl) % 2) + 1}_top{wildcards.percent}percent.bed",
        sample_bam_sorted="mm10_mapping/name_sorted/{AB}_{cond}_{repl}.bam",
    output:
        "peakqc/{AB}_{cond}_{repl,\d+}_top{percent}percent.csv"
    conda:
        '../env.yaml'
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
        echo "antibody,condition,replicate,topnpercent,mapped_fragments,mapped_fragments_in_peaks,frips,total_peaks,overlap_peaks,repl_percent" > {output}
        echo "{wildcards.AB},{wildcards.cond},{wildcards.repl},{wildcards.percent},$mappedFrags,$mappedFragsInPeaks,$frip,$total_peaks,$overlap_peaks,$repl_percent" >> {output}
    '''

rule plot_peak_qcs:
    input:
        csvs=expand('peakqc/{AB_cond}_{repl}_top{{percent}}percent.csv', AB_cond=AB_cond, repl=[1, 2]),
        beds=expand('peaks/{AB_cond}_{repl}_top{{percent}}percent.bed', AB_cond=AB_cond, repl=[1, 2])
    output:
        png='peakqc/plot_top{percent}percent.png',
        svg='peakqc/plot_top{percent}percent.svg',
    conda:
        '../env.yaml'
    script:
        '../scripts/plot_peak_qc.py'
