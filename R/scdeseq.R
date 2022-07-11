#' Run 'DESeq2' on aggregated data
#'
#' Run the differential expression analysis from 'DESeq2' on the aggregated
#' data. Users can filter out genes only present in a few samples.
#' @param agg_counts aggregated RNA counts
#' @param agg_table Table used to indicate which sample belongs to which state
#' @param design String used as formula in DE analysis (default = "~state")
#' @param comparison A string in form 'condition,group A,group B'
#' @param threads Number of workers used during DE analysis
#' @param do_prefilter Remove genes only expressed in a few samples
#' (default = TRUE)
#' @return Returns a DESeqDataSet object with results stored as metadata
#' columns.
#' @examples
#' # Basic run
#' example_dds <- scdeseq(agg_counts = agg_counts, agg_table = agg_table)
#'
#' # Do not prefilter
#' example_dds <- scdeseq(
#'   agg_counts = agg_counts, agg_table = agg_table,
#'   do_prefilter = FALSE
#' )
#'
#' @export
scdeseq <- function(agg_counts,
                    agg_table,
                    design = "~state",
                    comparison = "state,young,old",
                    threads = 10,
                    do_prefilter = TRUE) {

  # Define variables
  spl_com <- strsplit(comparison, split = ",")[[1]]
  con_id1 <- c(agg_table[, spl_com[1]] == spl_com[2])
  con_id2 <- c(agg_table[, spl_com[1]] == spl_com[3])
  design <- stats::as.formula(design)
  threads <- BiocParallel::register(BiocParallel::MulticoreParam(threads))

  # Generate DESeqDataSet for DE analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(agg_counts,
    colData = agg_table,
    design = design
  )

  # Remove the outliers
  if (do_prefilter) {
    dds_frag <- DESeq2::fpm(dds)

    kp_gene <- rowSums(dds_frag[, con_id1] >= 2) >= round(sum(con_id1) / 2) |
      rowSums(dds_frag[, con_id2] >= 2) >= round(sum(con_id2) / 2)

    de_out <- DESeq2::DESeq(dds[kp_gene, ], parallel = TRUE)
  } else {
    de_out <- DESeq2::DESeq(dds, parallel = TRUE)
  }

  # Return output
  return(de_out)
}
