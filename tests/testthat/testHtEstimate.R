library(crank) # We use the permute() function.

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
context("Dimensions")

probs = createProbMatrix(testmat1)

set.seed(4976401)
# Create sample outcome vector.
outcome = rnorm(testmat1_n)
# Choose an assignment vector for one of the permutations.
rand_column = sample(ncol(testmat1), 1)
assignment = testmat1[, rand_column]

# Compare assignment 1 to assignment 2.
result = htEstimate(outcome, assignment, c(1, -1, 0), probs)
