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

context("createProbMatrix - RI example 1")
y = c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
z = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
# Hide warning messages from the genperms function.
zz = capture.output({ perms = genperms(z, maxiter=10000) })
# This should not give us any errors.
prob_matrix = createProbMatrix(perms)
#prob_matrix

# TODO: what do we want to be checking here?
