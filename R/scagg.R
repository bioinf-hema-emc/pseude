#' Aggregate single-cell RNA-seq data to desired level
#'
#' This function allows the user to aggregate the counts from single-cell
#' RNA-seq experiments. Based on the metadata, cell can be linked to cell-types
#' (clusters), conditions (mutated vs wild type), etc. Aggregation will follow
#' these links.
#' @param sc_counts Table that contains the cell-level gene-counts
#' @param sc_table Look-up table that links cells to aggregation groups
#' @param sc_group String with column in agg_table that contains the grouping
#' @return Returns counts aggregated to specified level
#' @examples
#' # aggregate to sample level
#' example_agg <- scagg(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "sample_name"
#' )
#'
#' # aggregate to state level
#' example_agg <- scagg(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "state"
#' )
#'
#' @export
scagg <- function(sc_counts,
                  sc_table,
                  sc_group = "sample_name") {

  # alias external functions
  agg.Mat <- function(x, groupings = NULL, form = NULL, fun = "sum", ...) {
    Matrix.utils::aggregate.Matrix(x, groupings, form, fun, ...)
  }
  t.Mat <- function(x) {
    Matrix::t(x)
  }

  # aggregate the gene counts by summation; making a sample-level count-matrix
  return(t.Mat(agg.Mat(t.Mat(sc_counts), groupings = sc_table[, sc_group])))
}
