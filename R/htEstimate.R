#' @title Probabilities and joint probabilities of assignment.
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
      result[dest$start_row:dest$end_row, dest$start_col:dest$end_col] = t(targetMat)
    }
  }

  # Return the result.
  return(result)
}

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

#' @title Produce Horvitz-Thompson estimators of treatment assignment with standard error estimates, confidence intervals and hypothesis tests.
#'
#' @description TBD
#'
#' @param outcome Outcome vector for a given experiment.
#' @param assignment Assignment vector for the experiment.
#' @param contrasts A list of contasts. For example could be c(-1, 1, 0) for the above example if we wanted to compare treatment 1 to treatment 2. But in a factorial design, for example, we might want to compute the AMCE (see Hainmueller et al.). For example, if we had a 2⇥2⇥2 factorial design (8 treatment arms) and we wanted to look at an AMCE we might specify c(0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, - 0.25) and the elements of the treatment vector should be ti 2 {1,2,3,4,5,6,7,8}
#' @return estimate, standard_error, p value (two-tailed test of null)
#' @examples
#' TBD
htEstimate = function(outcome, raw_assignment, contrasts, prob_matrix, approx = "youngs", totals = F) {

  # Prepare some basic variables.

  # We sort the assignment levels in ascending order in case they are not consecutive natural numbers.
  assignment_levels = sort(unique(raw_assignment))

  # Create a dictionary mapping the raw assignment to the clean assignment levels.
  #assign_dict = list()
  #for (i in 1:length(assignment_levels)) {
  #  # assign_dict = assignment_levels[assignment_levels %in% assignment_levels]
  #  assign_dict[[assignment_levels[i]]] = i
 # }

  # Now create a clean_assignment field with consecutive natural numbers.
 # clean_assignment = assign_dict[[assignment]]

  # For now, just assume that all assignments are consecutive natural numbers (1:k).
  assignment = raw_assignment

  # K is the number of unique assignment values.
  k = length(assignment_levels)

  # N is the number of units.
  n = length(assignment)

  # 1. Calculate the effect estimate - see Aronow diss pages 14-15, esp. eq 2.8.

  # For each assignment, scale the outcomes of units with that assignment by the inverse of their
  # probability of assignment (on the diagonal matrix), generating a total for that outcome.

  outcome_totals = rep(NA, k)

  # Determine the location of the weights in the probability matrix.
  # We just need to identify the row because the column will be the same.
  weight_rows = c()
  for (i in 1:n) {
    weight_rows[i] = 1 + (i - 1) * n + (assignment[i] - 1)
  }

  # Loop over each assignment level.
  for (i in 1:k) {
    assignment_level = assignment_levels[i]
    weights = prob_matrix[weight_rows[assignment == assignment_level], weight_rows[assignment == assignment_level]]

    # Calculate the inverse-probability weighted total.
    outcome_totals[i] = outcome[assignment == assignment_level] / weights
  }

  # Then weight those outcomes by the contrast weightings, computing the sum, and divide by N.
  estimate = sum(outcome_totals * contrasts) / n

  # 2. Calculate the SE.

  # The sampling variance of the estimator is: (T = total, 1 = treatment, 0 = control)
  # 1/N^2 * (Var(Y1_T) + Var(Y0_T) - 2Cov(Y1_T, Y0_T))

  # Loop over each assignment level and calculate the variance of the total.
  # This is equation UEATE#21 (#27 when we account for clusters)
  # Actually, this is equation #32.
  variance_of_totals = rep(NA, k)
  covariances = matrix(nrow=k, ncol=k)
  for (assign_i in 1:k) {
    assignment_level = assignment_levels[assign_i]
    running_sum = 0

    for (i in 1:n) {
      # Pi_individual is the non-joint assignment probability of this unit to this level.
      # TODO: confirm that we should be using assign_i for the cells line, rather than just i.
      cells = getRawMatrixEntries(assign_i, assign_i, n)
      pi_i = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][i, i]

      # T'k(1 - pi_1k)*(Y_k/Pi_1k)^2
      # TODO: Does this need to check if pi_indiv == 0?? Doesn't seem to in the formula.
      first_component = (assignment[i] == assign_i) * (1 - pi_i)*(outcome[i]/pi_i)^2

      #individual_component = pi_individual * (1 - pi_individual) * (outcome[i]/pi_individual)^2
      running_sum = running_sum + first_component

      # Second and third components of equation 32
      for (j in 1:n) {
        # Skip this iteration if i == j
        if (i == j) next

        # TODO: confirm if we should use j or assign_i here.
        #cells = getRawMatrixEntries(assign_i, assign_i, n)
        pi_j = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][j, j]

        # TODO: confirm if we should use j or assign_i here.
        #cells = getRawMatrixEntries(assign_i, assign_i, n)
        pi_ij = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][i, j]

        if (pi_ij > 0) {
          # Option 1: pi_ij > 0
          # Second component of equation 32.
          joint_component = (assignment[i] == assign_i)*(assignment[j] == assign_i) / pi_ij *
            (pi_ij - pi_i * pi_j) * outcome[i] / pi_i * outcome[j] / pi_j
        } else {
          # Option 2: pi_ij = 0
          # Third component of equation 32.
          joint_component = (assignment[i] == assign_i) * outcome[i]^2/(2*pi_i)
            + (assignment[j] == assign_i) * outcome[j]^2/(2*pi_j)
        }

        running_sum = running_sum + joint_component
      }
    }
    variance_of_totals[assign_i] = running_sum

    # Create the covariances for this assignment level.
    # CUE#34
    cov_running_sum = 0
    # TODO: save work by only computing one triangle, rather than both triangles of the symmetric matrix.
    # TODO: right now this matrix is not symmetric, so there are remaining bugs.
    for (assign_j in 1:k) {
      # Skip when we are at our own level.
      if (assign_i == assign_j) next
      # Create a backup of the cells found in the previous loop, before we overwrite them.
      cells_i = cells
      for (i in 1:n) {

        # TODO: I think this part may be wrong, need to double-check. Should it be assign_i, assign_j?
        cells = getRawMatrixEntries(assign_i, assign_i, n)
        pi_i = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][i, i]
        for (j in 1:n) {
          # Skip the same observation.
          if (i == j) next

          # TODO: Should the cells indices be assign_i & assign_j? or i & j?
          cells = getRawMatrixEntries(assign_i, assign_j, n)
          pi_ij = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][i, j]

          cells = getRawMatrixEntries(assign_j, assign_j, n)
          pi_j = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][j, j]

          # TODO: need to review this equation, not sure if it's right.
          # Esp. the assignment[j] == assign_j part.
          # Equation #34 HERE (youngs inequality):
          # First component:
          cov_running_sum = cov_running_sum + (assignment[i] == assign_i) * (assignment[j] == assign_j) /
            pi_ij * (pi_ij - pi_i * pi_j) * outcome[i] * outcome[j] / (pi_i * pi_j)
        }

        # Eq#34 second component:
        cov_running_sum = cov_running_sum - (assign_i == assignment[i]) * outcome[i]^2 / (2 * pi_i)

        # Eq#34 third component:
        cells = getRawMatrixEntries(assign_j, assign_j, n)
        pi_0i = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][i, i]
        cov_running_sum = cov_running_sum - (assign_j == assignment[i]) * outcome[i]^2 / (2 * pi_0i)
      }
      covariances[assign_i, assign_j] = cov_running_sum

    }
  }

  # Weighted-sum of variance and covariance terms.
  # TODO: confirm that this use of the contrast weights is correct.
  var = sum(variance_of_totals * contrasts^2 + sum(contrasts %*% t(contrasts) * covariances, na.rm=T)) / n^2
  se = sqrt(var)

  # Calculate all needed covarianaces - we will use a weighted sum of their totals.

  # 3. Calculate the probability.

  # Return the results.
  result = list(estimate=estimate, std_err=se, p=NA)
  return(result)
}
