#' @title Probabilities and joint probailities of assignment.
#'
#' @description Return a matrix of probabilities and joint probabilities of assignment for n
#' observations and k randomization replications.
#'
#' @param assignments Each column is a replication and each element is the unit's assignment.
#' @param byrow Change the assignment matrix to use rows rather than columns for each randomization.
#' @return n*k by n*k matrix of probabilities and joint probabilities of assignment to each arm.
#' @examples
#' TBD
createProbMatrix = function(assignments, byrow = F) {
  # Return a n * k by n * k matrix.

  if (byrow) {
    # Transpose the assignment matrix if byrow is true.
    assignments = t(assignments)
  }

  # Determine the number of assignments by looking at the unique arms in the first replication.
  assign_levels = unique(assignments[, 1])
  k = length(assign_levels)
  n = nrow(assignments)

  # Matrix to store the results.
  result = matrix(nrow = n*k, ncol=n*k)
  #cat("Result dimensions:", n*k, "by", n*k, "\n")

  # Create the diagonals and lower triangle.
  # We start at the top-left diagonal and work down the column.
  for (col in 1:k) {
    for (row in col:k) {
      # Starting row & col in the matrix.
      start_row = 1 + (row - 1)  * n
      end_row = start_row + n - 1
      start_col = 1 + (col - 1) * n
      end_col = start_col + n - 1
      #cat("Start row:", start_row, "End row:", end_row, "Start col:", start_col, "End col:", end_col, "\n")
      probs = generateAssignmentProbs(row, col, assignments, assign_levels)
      #print(probs)
      result[start_row:end_row, start_col:end_col] = probs
    }
  }

  # Copy lower triangle to the upper triangle.

  # Return the result.
  return(result)
}

generateAssignmentProbs = function(row, col, assignments, assign_levels) {
  k = length(assign_levels)
  n = nrow(assignments)
  replications = ncol(assignments)

  # Our result will be an nxn matrix.
  probs = matrix(nrow = n, ncol = n)

  # Loop over the matrix.
  for (i in 1:n) {
    for (j in 1:n) {
      probs[j, i] = sum(assignments[i, ] == row & assignments[j, ] == col) / replications
    }
  }
  return(probs)
}

htEstimate = function(outcome, assignment, contrasts, prob_matrix, approx = "youngs", totals = F) {

}
