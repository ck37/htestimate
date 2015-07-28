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


# How does this compare to the ri package results?
library(ri)
#probs = genprobexact(assignment)
# These are not working, because RI only supports 0/1 assignment vectors :(
if (F) {
  testmat1
  probs = genprob(t(testmat1))
  probs
  probs = genprob(testmat1)
  probs
  perms = genperms(arms1)
  ate = estate(outcome, assignment, prob=probs)
  ate
}

# RI package example, but without blocking or clustering.
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
table(Z)
#cluster <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
#block <- c(rep(1,4),rep(2,6),rep(3,8))
#block

# Generates 10,000 samples by default.
perms = genperms(Z, maxiter=100000)
dim(perms)
# It stops at 43,758 because that is all of the permutations.
choose(18, 10)

# Look at the first 10 random permutations
perms[,1:10]
# Look at the assignment probability for those first 10 permutations.
colMeans(perms)[1:10]

probs <- genprob(perms) # probability of treatment
probs
ate <- estate(y,Z,prob=probs) # estimate the ATE
ate

# Convert from 0/1 to 1/2 assignment indictators, for compatability.
z_ck = Z + 1
table(z_ck)
perms_ck = perms + 1
dim(perms_ck)
# Double-check the first 10 permutations.
perms_ck[, 1:10]
prob_matrix = createProbMatrix(perms_ck)
# Examine the 1x1 matrix in the upper left.
prob_matrix[1:18, 1:18]

# Examine the htEstimate results.
result = htEstimate(y, z_ck, contrasts = c(-1, 1), prob_matrix = prob_matrix)
result

# Try with 100k perms
perms2 = genperms(Z, maxiter=100000)
perms_ck = perms2 + 1
prob_matrix2 = createProbMatrix(perms_ck)
result = htEstimate(y, z_ck, contrasts = c(-1, 1), prob_matrix = prob_matrix)
result


# TODO: Compare to the RI package example, this time with clustering and blocking:
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
cluster <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
block <- c(rep(1,4),rep(2,6),rep(3,8))
block
probs <- genprobexact(Z,blockvar=block, clustvar=cluster) # probability of treatment
probs
ate <- estate(y,Z,prob=probs) # estimate the ATE
ate
