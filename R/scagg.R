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
#' # Examples given with all optional parameters for education purposes
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
#' @export
scagg <- function(sc_counts,
                  sc_table,
                  sc_group = "sample_name") {

  # Deflate scounts ---------------------------
  if (class(sc_counts) != "dgCMatrix") {
    sc_counts <- Matrix::Matrix(as.matrix(sc_counts), sparse = TRUE)
  }

  # Transform grouping(s) and set formula ---------------------------
  A <- as.factor(do.call(paste, c(sc_table[, sc_group, drop = FALSE], sep = "_")))
  agg_group <- as.data.frame(A)

  agg_form <- stats::as.formula("~ 0 + A")

  # Generate model matrix for mapping ---------------------------
  agg_map <- Matrix::sparse.model.matrix(agg_form, agg_group, row.names = FALSE)

  colnames(agg_map) <- sub("^A", "", colnames(agg_map))

  # Aggregate by matrix multiplication ---------------------------
  return(sc_counts %*% agg_map)
}
