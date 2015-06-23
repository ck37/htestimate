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

# Check that the returned dimensions are correct.
context("Dimensions")

test_that("Matrix should have n * K rows", {
  result = createProbMatrix(testmat1)
  expect_equal(dim(result)[1], testmat1_n * testmat1_k)
})

# Matrix should have n x p columns.
test_that("Matrix should have n * K columns", {
  result = createProbMatrix(testmat1)
  expect_equal(dim(result)[2], testmat1_n * testmat1_k)
})


context("Diagonal")
# Check that the diagonal results are correct.

# Diagonal elements should sum to 1.

# Diagonal results should be symmetric.

# Diagonal results should vanish off-diagonal.

context("Off-Diagonal")
# Check that the off-diagonal results are correct.

# Check lower triangle.

# Confirm that base matrix is symmetric.
