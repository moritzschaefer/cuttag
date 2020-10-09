#!/bin/bash
# inspired by https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye.html

CHROM_SIZES=mm10.chrom_sizes

# mouse mapping
DNA-mapping --local --configFile ../config.yaml --input-dir fastq --output-dir mm10_mapping GRCm38_release98

# ecoli spikein mapping
DNA-mapping --local --configFile ../config.yaml --input-dir fastq --output-dir ecoli_mapping --alignerOpts '--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10' ecoli_K12

snakemake --use-conda -j 30 -s ../Snakefile


# mkdir normalized 2> /dev/null
# mkdir peaks 2> /dev/null

# SAMPLES=""
# for spikein_bam in ecoli_mapping/Bowtie2/*.bam; do
#     filename=$(basename -- "$spikein_bam")
#     mm10_bam=mm10_mapping/Bowtie2/$filename
#     sample=${filename%.bam}
#     SAMPLES="$sample $SAMPLES"

#     # normalize
#     seqDepthDouble=$(samtools view -F 0x04 "$spikein_bam" | wc -l)
#     seqDepth=$((seqDepthDouble/2))
#     if [[ "$seqDepth" -gt "1" ]]; then  # only continue if we have spike_in reads
#         scale_factor=`echo "10000 / $seqDepth" | bc -l`
#         echo "Scaling factor for $sample is: $scale_factor!"
#         # no need to filter for unmapped reads (they are not counted here)
#         bedtools genomecov -pc -bg -scale $scale_factor -ibam "$mm10_bam" > "normalized/$sample.bedgraph"

#     else
#         echo "seqDepth=0 for $sample"
#     fi
# done

# # seacr
# for sample in $SAMPLES

#         # call peaks, TODO: could also use relaxed instead of stringent
#         SEACR_1.3.sh "normalized/$sample.bedgraph" igg_control non stringent "peaks/${sample}_all"
#         SEACR_1.3.sh "normalized/$sample.bedgraph" 0.01 non stringent "peaks/${sample}_top1percent"



# # TODO get 1% data to review


# # SAMPLES='IgG_WT_2 IgG_Ago2_KO IgG_Ago21_KO IgG_Ago1_KO H3K9me3_WT_2 H3K9me3_Ago2_KO H3K9me3_Ago21_KO H3K9me3_Ago1_KO'
# mkdir mm10_mapping/sorted_bam 2>/dev/null
# echo "sample,repl,num_mapped_frags,num_mapped_frags_in_peaks,frips,num_peaks,num_repl_overlapped,repl_reproduced_percentage" > peak_qc.csv
# for sample in $SAMPLES; do
#     # calculate peak overlap across replicates
#     sample_short=${sample%_?}
#     repl=${sample##*_}
#     other_repl=$(( ( $repl % 2 ) + 1 ))

#     overlap_peaks=$(bedtools intersect -wa -a peaks/${sample}_all.stringent.bed -b peaks/${sample_short}_${other_repl}_all.stringent.bed | wc -l)
#     total_peaks=$(cat peaks/${sample}_all.stringent.bed | wc -l)
#     repl_percent=$(echo "100 * $overlap_peaks / $total_peaks" | bc -l)

#     # calculate FRagment proportion in Peaks regions (FRiPs)
#     mappedFragsDouble=$(samtools view -F 0x04 "mm10_mapping/Bowtie2/$sample.bam" | wc -l)
#     mappedFrags=$((mappedFragsDouble/2))
#     # weirdly, 'bedtools intersect' reports exactly twice as many reads as pairtobed. seems like bedtools intersect considers the whole fragment and always reports both reads per fragment
#     samtools sort -n "mm10_mapping/Bowtie2/$sample.bam" > mm10_mapping/sorted_bam/$sample.bam
#     mappedFragsInPeaks=$(bedtools pairtobed -type ospan -abam mm10_mapping/sorted_bam/$sample.bam -b peaks/${sample}_all.stringent.bed -bedpe | wc -l)
#     frip=$(echo "100 * $mappedFragsInPeaks / $mappedFrags" | bc -l)

#     echo "$sample,$repl,$mappedFrags,$mappedFragsInPeaks,$frip,$total_peaks,$overlap_peaks,$repl_percent" >> peak_qc.csv
# done

# # plot values
# python plot_peak_qc.py
