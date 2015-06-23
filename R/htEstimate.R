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

  #TODO: throw an error if k == 1
  #TODO: throw an error if n == 1

  # Matrix to store the results.
  result = matrix(nrow = n*k, ncol=n*k)
  #cat("Result dimensions:", n*k, "by", n*k, "\n")

  # Create the diagonals and lower triangle.
  # We start at the top-left diagonal and work down the column.
  for (col in 1:k) {
    for (row in col:k) {
      dest = getRawMatrixEntries(row, col, n)
      #cat("Start row:", start_row, "End row:", end_row, "Start col:", start_col, "End col:", end_col, "\n")
      probs = generateAssignmentProbs(row, col, assignments, assign_levels)
      #print(probs)
      result[dest$start_row:dest$end_row, dest$start_col:dest$end_col] = probs
    }
  }

  #cat("Current matrix after creating lower triangle:\n")
  #print(result)

  # Copy lower triangle to the upper triangle.

  for (col in 2:k) {
    for (row in seq(1, col-1)) {

      # Target is the matrix in the lower triangle that we're copying from.
      target = getRawMatrixEntries(col, row, n)
      #cat("Copying target:\n")
      #with(target, cat("Row:", col, "Col:", row, "Start row:", start_row, "End row:", end_row, "Start col:", start_col, "End col:", end_col, "\n"))

      targetMat = result[target$start_row:target$end_row, target$start_col:target$end_col]
      #print(targetMat)

      # Destination is the matrix in the upper triangle.
      dest = getRawMatrixEntries(row, col, n)

      # TODO: determine if we should be transposing the targetMat here.
      result[dest$start_row:dest$end_row, dest$start_col:dest$end_col] = targetMat
    }
  }

  # Return the result.
  return(result)
}

getRawMatrixEntries = function(row, col, n) {
  # Starting row & col in the matrix.
  start_row = 1 + (row - 1)  * n
  end_row = start_row + n - 1
  start_col = 1 + (col - 1) * n
  end_col = start_col + n - 1
  return(list(start_row=start_row, end_row=end_row, start_col=start_col, end_col=end_col))
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
