library(testthat) # For the context() function, in case we run this file manually.
library(crank) # We use the permute() function.
library(ri)

# Create test matrix 1 per the PDF document.
arms1 = 1:3
# We transpose the results so that each permutation is a column.
testmat1 = t(permute(arms1))
testmat1
# Number of units in the study.
testmat1_n = length(arms1)
# Number of assignments/arms.
testmat1_k = length(unique(arms1))

# Check that the returned dimensions are correct.
context("createProbMatrix - Dimensions")

test_that("Matrix should have n * K rows", {
  result = createProbMatrix(testmat1)
  expect_equal(dim(result)[1], testmat1_n * testmat1_k)
})

# Matrix should have n x p columns.
test_that("Matrix should have n * K columns", {
  result = createProbMatrix(testmat1)
  expect_equal(dim(result)[2], testmat1_n * testmat1_k)
})


# Confirm that there are no negative entries.
test_that("Matrix should have no negative entries", {
  result = createProbMatrix(testmat1)
  expect_equal(sum(apply(result, MARGIN=2, FUN=function(col) { sum(col < 0) } )), 0)
})

context("createProbMatrix - Diagonal")
# Check that the diagonal results are correct.

# Diagonal squares should sum to 1.

# Diagonal results should be symmetric.

# Diagonal results should vanish off-diagonal.

context("createProbMatrix - Off-Diagonal")
# Check that the off-diagonal results are correct.

# Each off-diagonal square should sum to 1 (I think).

# Check lower triangle.

# Confirm that base matrix is symmetric.

# Clean up test.
rm(arms1, testmat1, testmat1_n, testmat1_k)

context("createProbMatrix - RI example 1")
y = c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
z = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
# Hide warning messages from the genperms function.
zz = capture.output({ perms = genperms(z, maxiter=10000) })
# This should not give us any errors.
prob_matrix = createProbMatrix(perms)
#prob_matrix

# Confirm that there are no negative entries.
test_that("Matrix should have no negative entries", {
  expect_equal(sum(apply(prob_matrix, MARGIN=2, FUN=function(col) { sum(col < 0) } )), 0)
})

# TODO: what else do we want to be checking here?

# Clean up test.
rm(y, z, zz, prob_matrix, perms)

context("createProbMatrix - RI example 1")

context("createProbMatrix - CUE table 1")
# CUE example (table 1, p. 147), WITH blocking AND clustering.

ex = test_example_cue_table1()

# Review the Pi_1i's by cluster - probability of being assigned treatment.
# These should be 0.5 for the first 4 clusters then 0.333 for the last 6.
apply(ex$cluster_perms, MARGIN=1, FUN=function(x){ mean(x == "treatment") })

# We expect and do get 90 unique possible assignments.
dim(ex$cluster_perms)

# We should have a 16x90 matrix of assignment permutations at the unit level.
dim(ex$assign_perms)

# Confirm that our treatment assignment probabilities remain 0.5 and then 1/3.
apply(ex$assign_perms, MARGIN=1, FUN=function(x){ mean(x == "treatment") })
# We need to use this version once we have converted to numeric rather than string values.
apply(ex$assign_perms, MARGIN=1, FUN=function(x){ mean(x == 2) })

# This should give us a 32x32 matrix of assignment probabilities.
prob_matrix = createProbMatrix(ex$assign_perms)
dim(prob_matrix)

# Confirm that there are no negative entries.
test_that("Matrix should have no negative entries", {
  expect_equal(sum(apply(prob_matrix, MARGIN=2, FUN=function(col) { sum(col < 0) } )), 0)
})

# Clean up test.
rm(prob_matrix, ex)
