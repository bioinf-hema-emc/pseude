#' Aggregate single-cell RNA-seq data to desired level
#'
#' This function allows the user to aggregate the counts from single-cell
#' RNA-seq experiments. Based on the metadata, cell can be linked to cell-types
#' (clusters), conditions (mutated vs wild type), etc. Aggregation will follow
#' these links.
#' @param sc_counts Table that contains the cell-level gene-counts
#' @param sc_table Look-up table that links cells to aggregation groups
#' @param sc_group String with column in agg_table that contains the grouping
#' @param comparison A string in form "condition,group A,group B"
#' @return Returns a list with aggregated counts and aggregated annotation table to specified level
#' @examples
#' # Examples given with all optional parameters for education purposes
#' # aggregate to sample level
#' example_agg <- scagg(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "sample_name",
#'   comparison = "state,young,old"
#' )
#'
#' # aggregate to state level
#' example_agg <- scagg(
#'   sc_counts = raw_counts,
#'   sc_table = raw_table,
#'   sc_group = "state",
#'   comparison = "state,young,old"
#' )
#' @import Matrix
#' @export
scagg <- function(sc_counts,
                  sc_table,
                  sc_group = "sample_name",
                  comparison = "state,young,old") {

  # Deflate sc_counts ---------------------------
  if (class(sc_counts) != "dgCMatrix") {
    sc_counts <- Matrix::Matrix(as.matrix(sc_counts), sparse = TRUE)
  }

  # Subset sc_counts and sc_table ---------------------------
  spl_com <- strsplit(comparison, split = ",")[[1]]
  sc_table <- sc_table[is.element(sc_table[, spl_com[1]], spl_com[2:3]),]
  sc_counts <- sc_counts[, colnames(sc_counts) %in% rownames(sc_table)]

  # set comparison column to factor ---------------------------
  if (class(sc_table[, spl_com[1]]) != "factor") {
    sc_table[, spl_com[1]] <- as.factor(sc_table[, spl_com[1]])
    sc_table[, spl_com[1]] <- droplevels(sc_table[, spl_com[1]])
  }

  # Generate agg_table ---------------------------
  uni_table <- unique(sc_table[, c(sc_group, spl_com[1])])
  uni_fact <- as.factor(do.call(paste, c(uni_table[sc_group], sep = "_")))
  agg_table <- data.frame(uni_table[spl_com[1]], row.names = uni_fact)


  # Transform grouping(s) and set formula ---------------------------
  A <- as.factor(do.call(paste, c(sc_table[, sc_group, drop = FALSE], sep = "_")))
  agg_group <- as.data.frame(A)

  agg_form <- stats::as.formula("~ 0 + A")

  # Generate model matrix for mapping ---------------------------
  agg_map <- Matrix::sparse.model.matrix(agg_form, agg_group, row.names = FALSE)

  colnames(agg_map) <- sub("^A", "", colnames(agg_map))

  # Aggregate by matrix multiplication ---------------------------
  agg_list <- list(agg_counts = sc_counts %*% agg_map,
                   agg_table = agg_table)
  return(agg_list)
}
