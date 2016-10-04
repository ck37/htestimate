
#' @title TBD
#'
#' @description Used by createProbMatrix.
#'
#' @param row The assignment level to be used in the row check.
#' @param col The assignment level to be used in the column check.
#' @param assignments These have already been mapped to 1..k
#' @return nxn matrix with joint assignment probabilities.
generateAssignmentProbs = function(row, col, assignments) {
  n = nrow(assignments)
  replications = ncol(assignments)

  # Our result will be an nxn matrix.
  probs = matrix(nrow = n, ncol = n)

  # Loop over the matrix.
  for (i in 1:n) {
    for (j in 1:n) {
      # % of replications where Ai == row and Aj == col.
      probs[j, i] = sum(assignments[i, ] == row & assignments[j, ] == col) / replications
    }
  }
  return(probs)
}
