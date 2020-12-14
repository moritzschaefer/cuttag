# from Step28 on (https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye.html)
library(DESeq2) ## For differential analysis section
library(GenomicRanges)
library(dplyr)
library(chromVAR)

# Create a master peak list merging all the peaks called for each sample.
masterPeak = GRanges()
## overlap with bam file to get count
for(peaks in snakemake@input[["peaks"]]) {  # includes replicates
    peakRes = read.table(peaks, header = FALSE, fill = TRUE)
    masterPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(masterPeak, .)
}
masterPeak = reduce(masterPeak)

# save masterPeak file as bed
masterPeakDf <- data.frame(seqnames=seqnames(masterPeak),
                 starts=start(masterPeak)-1,
                 ends=end(masterPeak),
                 names=c(rep(".", length(masterPeak))),
                 scores=c(rep(".", length(masterPeak))),
                 strands=strand(masterPeak))
write.table(masterPeakDf, file=snakemake@output[["master_peaks"]], quote=F, sep="\t", row.names=F, col.names=F)


# Get the fragment counts for each peak in the master peak list.
library(DESeq2)
countMat = matrix(NA, length(masterPeak), length(snakemake@input[["bams"]]))
## overlap with bam file to get count
i = 1
for (bam in snakemake@input[["bams"]]) {
    fragment_counts <- getCounts(bam, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
}
colnames(countMat) = snakemake@params[["sample_names"]]

# Sequencing depth normalization and differential enriched peaks detection
selectR = which(rowSums(countMat) > 5) ## remove low count genes
dataS = countMat[selectR,]
condition = snakemake@params[["condition_names"]]
batch = snakemake@params[["batch_numbers"]]

colData = DataFrame(condition, batch)

save.image()

dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = colData,
                             design = ~batch + condition)  # we control for batch effect
dds = DESeq(dds)


saveRDS(dds, file=snakemake@output[["rds"]])
