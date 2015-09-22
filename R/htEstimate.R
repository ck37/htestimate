library(dplyr) # Used in htEstimate for the cluster aggregation.

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

  # Matrix to store the results. Initialize all cells to -1 so that we can see errors for easily.
  result = matrix(data=-1, nrow = n*k, ncol=n*k)
  #cat("Result dimensions:", n*k, "by", n*k, "\n")

  # Create the diagonals and lower triangle.
  # We start at the top-left diagonal and work down the column.
  for (col in 1:k) {
    for (row in col:k) {
      dest = getRawMatrixEntries(row, col, n)
      #cat("Start row:", start_row, "End row:", end_row, "Start col:", start_col, "End col:", end_col, "\n")
      probs = generateAssignmentProbs(row, col, assignments)
      #print(probs)
      result[dest$start_row:dest$end_row, dest$start_col:dest$end_col] = probs
    }
  }

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

#' @title  Unbiased treatment effect estimation with Horvitz-Thompson
#'
#' @description Produce Horvitz-Thompson estimators of treatment assignment with standard error estimates, confidence intervals
#' and hypothesis tests.
#'
#' @param outcome Outcome vector for a given experiment.
#' @param assignment Assignment vector for the experiment.
#' @param contrasts A list of contasts. For example could be c(-1, 1, 0) for the above example if we wanted to compare treatment 1 to treatment 2. But in a factorial design, for example, we might want to compute the AMCE (see Hainmueller et al.). For example, if we had a 2x2x2 factorial design (8 treatment arms) and we wanted to look at an AMCE we might specify c(0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, - 0.25) and the elements of the treatment vector should be ti 2 {1,2,3,4,5,6,7,8}
#' @param prob_matrix Probabilty matrix of assignment, as generated by createProbMatrix.
#' @param approx Options for bounding unidentified part of equation.
#' Default: "youngs" - Youngs inequality (see Aronow and Middleton).
#' Other options:
#' "constant effects" - constant effects assumption (See Aronow dissertation section 2.5),
#' "sharp null" - sharp null hypothesis (a special case of the constant effects assumption).
#' @param totals Calculate outcome totals rather than means, defaults to False.
#' @param cluster_id Cluster identifier, if data is in clusters. Outcomes will then be converted to cluster-totals.
#' @return estimate, standard_error, p value (two-tailed test of null)
#'
#' @examples
#'
#' # Example using data from RI package.
#' y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
#' z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
#' # Generate 10,000 random permutations of the assignment vector.
#' perms = ri::genperms(z, maxiter=10000)
#' # Estimate the probability of assignment for each unit and assignment level.
#' prob_matrix = createProbMatrix(perms)
#' # Estimate the treatment effect using Horvitz-Thompson.
#' htestimate(y, z, contrasts = c(-1, 1), prob_matrix = prob_matrix)
#'
htestimate = function(outcome, assignment, contrasts, prob_matrix, approx = "youngs", totals = F, cluster_id = NULL) {

  # Rename the original version of the assignments before we renumber the assignment levels (potentially).
  raw_assignment = assignment

  # Prepare some basic variables.

  # Handle clustering if the cluster id was defined.
  if (!is.null(cluster_id)) {
    cat("Handling clustering with cluster_id variable.\n")

    # Create a combined dataframe for use with dplyr.
    combined_df = data.frame(outcome, raw_assignment, cluster_id)

    # Aggregate the data by cluster id.
    agg_df = dplyr::group_by(combined_df, cluster_id)

    # Confirm that all assignments are equal within each cluster, otherwise give an error.
    assignments_by_cluster = dplyr::n_distinct(agg_df, raw_assignment)
    count_of_nonunique_clusters = sum(assignments_by_cluster == 1)
    if (count_of_nonunique_clusters != 0) {
      throw("Assignments differ within clusters.")
    }

    # TODO: Aggregate the outcome data into cluster totals.

    # TODO: Aggregate the assignment data into cluster totals.

    # TODO: do we need to mess with the probability matrix? Collapse it to clusters?

  }

  # We sort the assignment levels in ascending order in case they are not consecutive natural numbers.
  raw_assignment_levels = sort(unique(raw_assignment))

  # Create a dictionary mapping the raw assignment to the clean assignment levels.
  assign_dict = list()
  for (i in 1:length(raw_assignment_levels)) {
    # assign_dict = assignment_levels[assignment_levels %in% assignment_levels]
    #cat("Assign levels", i, "has", assign_levels[[i]])
    # Convert it to a character because 0 cannot be used as a list name in R.
    assign_dict[[as.character(raw_assignment_levels[i])]] = i
  }

  # Now create a clean assignment field with consecutive natural numbers.
  # Here we just need to loop over elements, because it's an array rather than a matrix as in createProbMatrix.
  # That's why we use sapply here.
  assignment = sapply(raw_assignment, FUN=function(x) { assign_dict[[as.character(x)]] })

  # NOTE: createProbMatrix does a similar operation to convert assignment levels to natural numbers.

  # Now use the new assignment levels within the function.
  assignment_levels = unlist(assign_dict, use.names=F)

  # K is the number of unique assignment values.
  k = length(assignment_levels)

  # N is the number of units.
  n = length(assignment)
  #cat("N:", n, "K:", k, "\n")

  # 1. Calculate the effect estimate - see Aronow diss pages 14-15, esp. eq 2.8.

  # For each assignment, scale the outcomes of units with that assignment by the inverse of their
  # probability of assignment (on the diagonal matrix), generating a total for that outcome.

  outcome_totals = rep(NA, k)

  # Determine the location of the weights in the probability matrix.
  # We just need to identify the row because the column will be the same.
  weight_rows = c()
  for (i in 1:n) {
    weight_rows[i] = 1 + (i - 1) * n + (assignment[i] - 1)
  #  cat("i:", i, "assignment[i]", assignment[i], "Weight index:", weight_rows[i], "\n")
  }
  #cat("Weight_rows:", weight_rows, "\n")

  # Loop over each assignment level.
  #cat("Loop over each assignment level.\n")
  for (i in 1:k) {
    assignment_level = assignment_levels[i]
    # TODO: are these the right cells in the prob matrix? Doesn't it need to vary by the assignment level?
    # I think we need to use getRawMatrixEntries function here, based on the assignment level.
    #weights = prob_matrix[weight_rows[assignment == assignment_level], weight_rows[assignment == assignment_level]]
    cells = getRawMatrixEntries(assignment_level, assignment_level, n)
    #cat("i", i, "assignment_level", assignment_level, "Cells:", "\n")
    #print(cells)
    weights = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][assignment == assignment_level, assignment == assignment_level]
    # Number of observations with this assignment level
    n_obs = length(outcome[assignment == assignment_level])
    #cat("Weights:\n")
    #print(weights)
    #cat("Dim weights:\n")
    #print(dim(weights))
    #cat("Outcomes: (", n_obs, ")\n")
    #print(outcome[assignment == assignment_level])
    #cat("Diag of weights:\n")
    #print(matDiag(weights))

    # Calculate the inverse-probability weighted total.
    # This doesn't work correctly - commenting out and running in a loop for now.
    # outcome_totals[i] = outcome[assignment == assignment_level] / matDiag(weights)
    # Do it slower for debugging.
    # TODO: convert to dot product?
    outcome_totals[i] = 0
    for (diag_i in 1:n_obs) {
      outcome_totals[i] = outcome_totals[i] + outcome[assignment == assignment_level][diag_i] / matDiag(weights)[diag_i]
    }
    #cat("Outcome totals:\n")
    #print(outcome_totals[i])
  }

  #cat("Completed the loop.\n")
  #cat("Outcome totals:", outcome_totals, "\n")

  # Then weight those outcomes by the contrast weightings, computing the sum, and divide by N.
  if (totals) {
    # Calculate the difference in totals.
    # TODO: confirm this is the correct way with Joel.
    estimate = sum(outcome_totals * contrasts)
  }
  else {
    # Calculate the difference in means (default).
    estimate = sum(outcome_totals * contrasts) / n
  }
  #cat("Estimate:", estimate, "\n")

  # 2. Calculate the SE.

  # Calculate all needed covarianaces - we will use a weighted sum of their totals.


  # Loop over each assignment level and calculate the variance of the total.
  # Actually, this is equation CUEATE#32.
  # TODO: generalize EQ#32 to multiple treatment arms for reference.
  variance_of_totals = rep(NA, k)
  covariances = matrix(nrow=k, ncol=k)
  for (assign_a in 1:k) {
    # TODO: use level_a when the assignment levels are not 1:k
    # Actually - we have already renumbered the assignment levels so this is probably not needed.
    level_a = assignment_levels[assign_a]
    running_sum = 0
    #cat("Assign_a:", assign_a, "Level_a:", level_a, "\n")

    # Loop over each observation.
    for (obs_k in 1:n) {
      #cat("Assignment level:", level_a, "Observation:", obs_k)
      # Pi_ai is the non-joint assignment probability of this unit to this assignment level/arm.
      # TODO: confirm that we should be using assign_i for the cells line, rather than just i.
      # TODO: should we be using level_a here?
      cells = getRawMatrixEntries(assign_a, assign_a, n)
      #cat("Cells:\n")
      #print(cells)
      pi_ak = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_k, obs_k]

      # Equation: UEATE#32, first component.
      # T'k(1 - pi_1k)*(Y_k/Pi_1k)^2
      # TODO: Does this need to check if pi_indiv == 0?? Doesn't seem to in the formula.
      first_component = (assignment[obs_k] == assign_a) * (1 - pi_ak)*(outcome[k]/pi_ak)^2

      #individual_component = pi_individual * (1 - pi_individual) * (outcome[i]/pi_individual)^2
      running_sum = running_sum + first_component

      # Second and third components of equation 32
      for (obs_l in 1:n) {
        # Skip this iteration if i == j
        if (obs_k == obs_l) next

        # TODO: confirm if we should use j or assign_i here.
        #cells = getRawMatrixEntries(assign_i, assign_i, n)
        pi_al = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_l, obs_l]

        # TODO: confirm if we should use j or assign_i here.
        #cells = getRawMatrixEntries(assign_i, assign_i, n)
        pi_ak_al = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_k, obs_l]

        if (pi_ak_al > 0) {
          # Option 1: pi_ij > 0
          # Second component of equation 32.
          joint_component = (assignment[obs_k] == assign_a)*(assignment[obs_l] == assign_a) / pi_ak_al *
            (pi_ak_al - pi_ak * pi_al) * outcome[obs_k] / pi_ak * outcome[obs_l] / pi_al
        } else {
          # Option 2: pi_ij = 0
          # Third component of equation 32.
          # This is young's equality right now.
          # TODO: support constant effects variance estimation (Aronow diss 2.5)
          # TODO: support sharp null hypothesis.
          if (approx == "youngs") {
            joint_component = (assignment[obs_k] == assign_a) * outcome[obs_k]^2/(2*pi_ak)
              + (assignment[obs_l] == assign_a) * outcome[obs_l]^2/(2*pi_al)
          } else if (approx == "constant effects") {
            # TODO: implementation here.
          } else if (approx == "sharp null") {
            # TODO: implementation here.
          }
        }

        running_sum = running_sum + joint_component
      }
    }
    variance_of_totals[assign_a] = running_sum

    # Create the covariances for this assignment level.
    # CUE#34
    cov_running_sum = 0
    # TODO: save work by only computing one triangle, rather than both triangles of the symmetric matrix.
    # TODO: right now this matrix is not symmetric, so there are remaining bugs.
    for (assign_b in 1:k) {
      # Skip when we are at our own level.
      if (assign_a == assign_b) next

      level_b = assignment_levels[assign_b]
      # Create a backup of the cells found in the previous loop, before we overwrite them.
      cells_i = cells
      for (obs_k in 1:n) {

        # TODO: I think this part may be wrong, need to double-check. Should it be assign_i, assign_j?
        cells = getRawMatrixEntries(assign_a, assign_a, n)
        pi_ak = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_k, obs_k]
        for (obs_l in 1:n) {
          # Skip the same observation.
          if (obs_k == obs_l) next

          # TODO: Should the cells indices be assign_i & assign_j? or i & j?
          cells = getRawMatrixEntries(assign_a, assign_b, n)
          pi_ak_bl = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_k, obs_l]

          cells = getRawMatrixEntries(assign_b, assign_b, n)
          pi_bl = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_l, obs_l]

          # TODO: need to review this equation, not sure if it's right.
          # Esp. the assignment[j] == assign_j part.
          # Equation #34 HERE (youngs inequality):
          # First component:
          cov_running_sum = cov_running_sum + (assignment[obs_k] == assign_a) * (assignment[obs_l] == assign_b) /
            pi_ak_bl * (pi_ak_bl - pi_ak * pi_bl) * outcome[obs_k] * outcome[obs_l] / (pi_ak * pi_bl)
        }

        # Eq#34 second component:
        cov_running_sum = cov_running_sum - (assign_a == assignment[obs_k]) * outcome[obs_k]^2 / (2 * pi_ak)

        # Eq#34 third component:
        cells = getRawMatrixEntries(assign_b, assign_b, n)
        pi_bk = prob_matrix[cells$start_row:cells$end_row, cells$start_col:cells$end_col][obs_k, obs_k]
        cov_running_sum = cov_running_sum - (assign_b == assignment[obs_k]) * outcome[obs_k]^2 / (2 * pi_bk)
      }
      covariances[assign_a, assign_b] = cov_running_sum

    }
  }

  # The sampling variance of the estimator is: (T = total, 1 = treatment, 0 = control)
  # 1/N^2 * (Var(Y1_T) + Var(Y0_T) - 2Cov(Y1_T, Y0_T))

  # Weighted-sum of variance and covariance terms.
  # TODO: confirm that this use of the contrast weights is correct.
  # The general formula is from https://en.wikipedia.org/wiki/Variance#Weighted_sum_of_variables
  # I think we don't need to multiply two because we are using both triangles of a symmetric matrix.
  var = sum(variance_of_totals * contrasts^2 + sum(contrasts %*% t(contrasts) * covariances, na.rm=T))

  # Divide by n^2 if we're not calculating the VAR of totals.
  if (!totals) {
     var = var / n^2
  }
  se = sqrt(var)

  # 3. Calculate the probability using the normal distribution.
  # TODO: could we use RI or something else as another way to generate the p-value?
  p_value = 2*pnorm(-abs(estimate/se))

  # Return the results.
  # We include the variances and covariances for debugging purposes only.
  result = list(estimate=estimate, std_err=se, p_value=p_value, variances=variance_of_totals, covariances=covariances)
  return(result)
}

#' @title Return the diagonal of the matrix.
#' @description Returns the diagonal elements of a matrix as a vector.
#' @param mat Input matrix, must be square.
#' @return diagonal elements of the matrix as a vector.
#' @examples
#' myMatrix = matrix(1:9, ncol=3)
#' # This returns c(1, 5, 9)
#' diagonal = matDiag(myMatrix)
matDiag = function(mat) {
  # Must be a 1x1 vector.
  if (class(mat) != "matrix") {
    if (length(mat) == 1) {
      return(mat[1])
    } else {
      throw("Input argument is not a matrix")
    }
  }
  if (dim(mat)[1] != dim(mat)[2]) {
    throw("Diagonal is only defined on square matrices.")
  }

  # Allocate our results first so that we don't have to resize the vector.
  diag = rep(0, dim(mat)[1])

  for (i in 1:dim(mat)[1]) {
    diag[i] = mat[i, i]
  }
  return(diag)
}

# New function to randomly permute assignments within blocks, after aggregating to cluster totals.
# Or do randomizr, blockTools, or ri (et al?) already provide this functionality?
