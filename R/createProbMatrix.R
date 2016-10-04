#' @title Probabilities and joint probabilities of assignment.
#'
#' @description Return a matrix of probabilities and joint probabilities of assignment for n
#' observations and k randomization replications.
#'
#' @param assignments Each column is a replication and each element is the unit's assignment.
#' @param byrow Change the assignment matrix to use rows rather than columns for each randomization.
#' @return n*k by n*k matrix of probabilities and joint probabilities of assignment to each arm.
#' @examples
#' # Create test matrix 1 per the PDF document.
#' arms1 = 1:3
#' # We transpose the results so that each permutation is a column.
#' testmat1 = t(crank::permute(arms1))
#' # Display the simple test matrix.
#' testmat1
#' prob_matrix = createProbMatrix(testmat1)
#'
#' # RI package example, but without blocking or clustering.
#' y = c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
#' z = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
#' # Create 10k random permutations of the assignment vector.
#' perms = ri::genperms(z, maxiter=10000)
#' prob_matrix = createProbMatrix(perms)
#'
createProbMatrix = function(assignments, byrow = F) {
  # Return a n * k by n * k matrix.

  # Rename the original version of the assignments before we re-number the assignment levels (potentially).
  raw_assignments = assignments

  if (byrow) {
    # Transpose the assignment matrix if byrow is true.
    raw_assignments = t(raw_assignments)
  }

  # Convert assignment levels to be 1..k if they aren't already (e.g. if they are 0/1 per ri package).
  # Here we take the first column only to generate the assignment levels.
  assignment_levels = sort(unique(raw_assignments[,1]))

  # Create a dictionary mapping the raw assignment to the clean assignment levels.
  assign_dict = list()
  for (i in 1:length(assignment_levels)) {
    # assign_dict = assignment_levels[assignment_levels %in% assignment_levels]
    #cat("Assign levels", i, "has", assign_levels[[i]])
    # Convert it to a character because 0 cannot be used as a list name in R.
    assign_dict[[as.character(assignment_levels[i])]] = i
  }

  # Now create a clean assignment field with consecutive natural numbers.
  assignments = apply(raw_assignments, MARGIN=c(1, 2), FUN=function(x) { assign_dict[[as.character(x)]] })

  # Determine the number of assignments by looking at the unique arms in the first replication.
  #assign_levels = unique(assignments[, 1])
  k = length(assignment_levels)
  n = nrow(assignments)

  if (k == 1) {
    stop("Error: found only 1 assignment level; should be 2 or greater.")
    # TODO: add a test case to confirm that this check works.
  }

  if (n == 1) {
    stop("Error: only 1 observation found; should be 2 or greater.")
    # TODO: add a test case to confirm that this check works.
  }

  # Create the stacked indicator matrices.
  stacked_inds = matrix(nrow = n*k, ncol = ncol(assignments))

  for (assign in 1:k) {
    indicator_matrix = as.numeric(assignments == assign)
    stacked_inds[(n*(assign-1)+1):(n*assign), ] = indicator_matrix
  }

  # Use the stacked indicator matrices to calculate the probability matrix.
  result = stacked_inds %*% t(stacked_inds) / ncol(assignments)

  # Label the rows and columns based on the raw assignment labels (or codes if there are no labels).

  # Column name should be: [assign col]_[obs #]
  col_names = paste(rep(paste(assignment_levels), each=n), 1:n, sep="_")
  stopifnot(length(col_names) == ncol(result))
  colnames(result) = col_names
  # Row name should be: [assign row]_[obs #]
  row_names = paste(rep(paste(assignment_levels), each=n), 1:n, sep="_")
  stopifnot(length(row_names) == nrow(result))
  rownames(result) = row_names


  # Return the result.
  return(result)
}
