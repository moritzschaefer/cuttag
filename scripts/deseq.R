# from Step28 on (https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye.html)
library(DESeq2) ## For differential analysis section
library(GenomicRanges)
library(dplyr)
library(chromVAR)

# Create a master peak list merging all the peaks called for each sample.
mPeak = GRanges()
## overlap with bam file to get count
for(peaks in snakemake@input[["peaks"]]) {  # includes replicates
    peakRes = read.table(peaks, header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}
masterPeak = reduce(mPeak)

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

save.image()

dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

countMatDiff = cbind(dataS, normDDS, res)
head(countMatDiff)

# copied from rna-seq-star-deseq2/deseq2.R
res <- res[order(res$padj),]

svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]], sep=",")
