#' A wrapper around the four main functions
#'
#' This function allows the user to run through the entire pipeline in one go.
#' The pipeline contains scagg, scdeseq, scout and schow. This does reduces the
#' parameter-tuning; For specifics, please use the separate functions
#' @param scounts Table that contains the cell-level gene-counts
#' @param scable Look-up table that links cells to aggregation groups
#' @param scroup String with column in agg_table that contains the grouping
#' @param design String used as formula in DE analysis (default = "~state")
#' @param comparison A string in form 'condition,group A,group B'
#' @param threads Number of workers used during DE analysis
#' @param do_prefilter Remove genes only expressed in a few samples
#' @param retorwrite Binary choice, return object or write to files
#' @param out_path Which dir should house the output files
#' @return Optional: Returns DESeq2 data set with attached results and
#' information for visualization in a list
#' @examples
#' # Return object list
#' example_out <- scuard(
#'   scounts = raw_counts,
#'   scable = raw_table,
#'   scroup = "sample_name",
#'   design = "~state",
#'   comparison = "state,young,old",
#'   out_path = NULL,
#'   threads = 10,
#'   do_prefilter = TRUE,
#'   retorwrite = TRUE
#' )
#' # Write to files
#' scuard(
#'   scounts = raw_counts,
#'   scable = raw_table,
#'   scroup = "sample_name",
#'   design = "~state",
#'   comparison = "state,young,old",
#'   threads = 10,
#'   do_prefilter = TRUE,
#'   retorwrite = FALSE
#' )
#' @export
scuard <- function(scounts,
                   scable,
                   scroup = "sample_name",
                   design = "~state",
                   comparison = "state,young,old",
                   out_path = NULL,
                   threads = 10,
                   do_prefilter = TRUE,
                   retorwrite = TRUE) {

  # Call scagg
  aggounts <- scagg(scounts, scable, scroup)

  # Create aggable
  A <- as.factor(do.call(paste, c(test_table[, test_group, drop = FALSE], sep = "_")))


  # Call scdeseq
  aggdds <- scdeseq(aggounts, aggable, design,  comparison, threads, do_prefilter)


  }
