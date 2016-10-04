#' @title Return probability matrix locations for the sub-matrix for a certain observation pair.
#'
#' @description TBD
#'
#' @param row TBD
#' @param col TBD
#' @param n Total number of observations.
#' @param direct Not yet implemeneted: try to return the ranges for direct usage, rather than a list.
#' @return TBD
# Doesn't this need the number of assignments and/or the total number of records? -> No, just n.
getRawMatrixEntries = function(row, col, n, direct=F) {
  # Starting row & col in the full matrix.
  start_row = 1 + (row - 1)  * n
  end_row = start_row + n - 1
  start_col = 1 + (col - 1) * n
  end_col = start_col + n - 1
  if (direct) {
    # TODO: can we do something like this?
    #return(c(c(start_row:end_row), c(start_col:end_col)))
  } else {
    return(list(start_row=start_row, end_row=end_row, start_col=start_col, end_col=end_col))
  }
}
