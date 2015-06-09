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

  # Determine the number of assignments.
  assign_levels = unique(assignments)
  k = length(assign_levels)
  n = nrow(assignments)

  # First, create the diagonal elements.

  # Then create the upper triangle.

  # Copy upper triangle to the lower triangle.
}

htEstimate = function(outcome, assignment, contrasts, prob_matrix, approx = "youngs", totals = F) {

}
