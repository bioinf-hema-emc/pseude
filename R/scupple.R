#' A wrapper around the four main functions
#'
#' This function allows the user to run through the entire pipeline in one go.
#' The pipeline contains scagg, sciff, scout and scuality. This does reduces the
#' parameter-tuning; For specifics, please use the separate functions
#' @param sc_counts Table that contains the cell-level gene-counts
#' @param sc_table Look-up table that links cells to aggregation groups
#' @param sc_group String with column in agg_table that contains the grouping
#' @param design String used as formula in DE analysis (default = "state")
#' @param comparison A string in form 'condition,group A,group B'
#' @param out_path Which dir should house the output files
#' @param threads Number of workers used during DE analysis
#' @param do_prefilter Remove genes only expressed in a few samples
#' @param writeorret Binary choice, write to files or return object
#' @return Optional: Returns DESeq2 data set with attached results and
#' information for visualization in a list
#' @examples
#' # Examples given with all optional parameters for education purposes
#' \dontrun{
#' # Write to files
#' scupple(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "sample_name",
#'   design = "state",
#'   comparison = "state,young,old",
#'   out_path = "/example/output/path",
#'   threads = 10,
#'   do_prefilter = TRUE,
#'   writeorret = TRUE
#' )
#' }
#' # Return object list
#' example_out <- scupple(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "sample_name",
#'   design = "state",
#'   comparison = "state,young,old",
#'   out_path = NULL,
#'   threads = 10,
#'   do_prefilter = TRUE,
#'   writeorret = FALSE
#' )
#' @export
scupple <- function(sc_counts,
                   sc_table,
                   sc_group = "sample_name",
                   design = "state",
                   comparison = "state,young,old",
                   out_path = NULL,
                   threads = 10,
                   do_prefilter = TRUE,
                   writeorret = TRUE) {

  # Call scagg ---------------------------
  agg_list <- scagg(
    sc_counts = sc_counts,
    sc_table = sc_table,
    sc_group = sc_group,
    comparison = comparison
  )

  # Call scdeseq ---------------------------
  agg_dds <- sciff(
    agg_counts = agg_list$agg_counts,
    agg_table = agg_list$agg_table,
    design = design,
    comparison = comparison,
    threads = threads,
    do_prefilter = do_prefilter
  )

  # Test for writeorret, call scout and schow ---------------------------
  if (writeorret) {
    scout(
      dds = agg_dds,
      comparison = comparison,
      out_path = out_path,
      do_write = TRUE,
      do_return = FALSE
    )
    scuality(
      dds = agg_dds,
      comparison = comparison,
      out_path = out_path,
      do_write = TRUE,
      do_labels = TRUE,
      do_rlog = FALSE,
      do_return = FALSE
    )
  } else {
    agg_out <- scout(
      dds = agg_dds,
      comparison = comparison,
      out_path = out_path,
      do_write = FALSE,
      do_return = TRUE
    )
    agg_show <- scuality(
      dds = agg_dds,
      comparison = comparison,
      out_path = out_path,
      do_write = FALSE,
      do_labels = TRUE,
      do_rlog = FALSE,
      do_return = TRUE
    )
    return(list(agg_out, agg_show))
  }
}
