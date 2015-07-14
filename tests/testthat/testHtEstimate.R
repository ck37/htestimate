library(crank) # We use the permute() function.
library(testthat)

# Create test matrix 1 per the PDF document.
arms1 = 1:3
# We transpose the results so that each permutation is a column.
testmat1 = t(permute(arms1))
testmat1
# Number of units in the study.
testmat1_n = length(arms1)
# Number of assignments/arms.
testmat1_k = length(unique(arms1))

# Check that a simple test case works.
context("Simple test case")

prob_matrix = createProbMatrix(testmat1)
prob_matrix

set.seed(4976401)
# Create sample outcome vector.
outcome = rnorm(testmat1_n)
outcome

# Choose an assignment vector for one of the permutations.
rand_column = sample(ncol(testmat1), 1)
assignment = testmat1[, rand_column]
assignment
# Make a copy equal to the internal argument of htEstimate, to help with debugging.
raw_assignment = assignment

# Compare assignment 1 to assignment 2.
contrasts = c(1, -1, 0)
approx = "youngs"
result = htEstimate(outcome, assignment, contrasts, prob_matrix)
result

# This is what the variance weights are for the standard error estimate.
contrasts %*% t(contrasts)
