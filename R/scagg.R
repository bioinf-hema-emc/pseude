#' Aggregate single-cell RNA-seq data to desired level
#'
#' This function allows the user to aggregate the counts from single-cell
#' RNA-seq experiments. Based on the metadata, cell can be linked to cell-types
#' (clusters), conditions (mutated vs wild type), etc. Aggregation will follow
#' these links.
#' @param scounts Table that contains the cell-level gene-counts
#' @param scable Look-up table that links cells to aggregation groups
#' @param scroup String with column in agg_table that contains the grouping
#' @return Returns counts aggregated to specified level
#' @examples
#' # aggregate to sample level
#' example_agg <- scagg(
#'   scounts = raw_counts,
#'   scable = raw_table,
#'   scroup = "sample_name"
#' )
#'
#' # aggregate to state level
#' example_agg <- scagg(
#'   scounts = raw_counts,
#'   scable = raw_table,
#'   scroup = "state"
#' )
#' @export
scagg <- function(scounts,
                  scable,
                  scroup = "sample_name") {

  # Deflate sc_counts ---------------------------
  if (class(x) != "dgCMatrix") {
    scounts <- Matrix::Matrix(as.matrix(scounts), sparse = TRUE)
  }

  # Transform grouping(s) and set formula ---------------------------
  A <- as.factor(do.call(paste, c(scable[, scroup, drop = FALSE], sep = "_")))
  aggroup <- as.data.frame(A)

  aggform <- stats::as.formula("~ 0 + A")

  # Generate model matrix for mapping ---------------------------
  aggmap <- Matrix::sparse.model.matrix(aggform, aggroup, row.names = FALSE)

  colnames(aggmap) <- sub("^A", "", colnames(aggmap))

  # Aggregate by matrix multiplication ---------------------------
  return(scounts %*% aggmap)
}
