library("DESeq2")

dds <- readRDS(snakemake@input[["rds"]])
dataS = counts(dds, normalized = FALSE)
normDDS = counts(dds, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

res = results(dds, independentFiltering = FALSE,
              altHypothesis = "greaterAbs", alpha=0.05,
              contrast=c("condition", snakemake@wildcards[["target"]], snakemake@params[["control"]])
              )


countMatDiff = cbind(dataS, normDDS, res)

# copied from rna-seq-star-deseq2/deseq2.R
res <- res[order(res$padj),]

svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

rownames(countMatDiff) <- lapply(1:nrow(countMatDiff), function(i) paste0("masterPeak_", i))

write.table(countMatDiff, file=snakemake@output[["table"]], sep=",")
