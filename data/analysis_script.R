# Analysis script

library(pseude)



agg_counts <- scaggregate(raw_counts, agg_dt)

agg_tell <- 1

if (class(agg_tell) == "character") {
  annot_cols <- colnames(agg_dt)[colnames(agg_dt) != agg_tell]
} else {
  annot_cols <- colnames(agg_dt)[-agg_tell]
}

unique(agg_dt[order(sample_name), annot_cols, with = FALSE])



annot_cols



annot_table <- unique(agg_dt[order(sample_name)])


design <- as.formula("~state")

dds <- DESeqDataSetFromMatrix(agg_counts, colData = annotation_table, design = design)

de_out <- DESeq2::DESeq(dds, parallel = FALSE)
