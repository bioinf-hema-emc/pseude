#' Aggregate single-cell RNA-seq data to desired level
#' 
#' This function allows the user to aggregate the counts from single-cell 
#' RNA-seq experiments. Based on the metadata, cell can be linked to cell-types
#' (clusters), conditions (mutated vs wild type), etc. Aggregation will follow 
#' these links.
#' @param sc_counts Table that contains the cell-level gene-counts
#' @param agg_link Look-up table that links cells to aggregation groups
#' @param agg_group String with column in agg_link that contains the grouping
#' @return Returns counts aggregated to specified level
#' @examples 
#' example_mat <- matrix(sample(0:15, 500 * 240, TRUE), 500, 240)
#' colnames(example_mat) <- rownames(example_data)
#' rownames(example_mat) <- paste("gene", seq(1:500), sep = "_")
#' 
#' example_condition <- rep(rep(c("BM", "CTC"), each = 20), 6)
#' 
#' example_patient <- rep(c(
#'   "p2481", "p2482", "p5312", "p8945", "p8765",
#'   "p1264"
#' ), each = 20)
#' 
#' example_link <- data.frame(
#'   sample_ID = paste(example_patient,
#'     example_condition,
#'     sep = "_"
#'   ),
#'   condition = example_condition,
#'   patient_ID = example_patient,
#'   row.names = paste("cell", seq(1:240), sep = ".")
#' )
#' 
#' example_agg <- scaggregate(sc_counts = example_mat, agg_link = example_data, agg_group)
#' 
#' @export
scaggregate <- function(sc_counts, agg_link, agg_group = "sample_name") {
  
  # alias external functions
  agg.Mat <- function(x, groupings = NULL, form = NULL, fun = "sum", ...) {
    Matrix.utils::aggregate.Matrix(x, groupings, form, fun, ...)
  }
  t.Mat <- function(x) {
    Matrix::t(x)
  }
  
  # Check whether sc_counts is sparse
  if (class(sc_counts) == 'dgCMatrix') {
    # aggregate the gene counts by summation; making a sample-level count-matrix
    return(t.Mat(agg.Mat(t.Mat(sc_counts), groupings = agg_link[, agg_group])))
  } else {
    return(NULL)
  }
}