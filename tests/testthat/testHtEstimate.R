library(crank) # We use the permute() function.
library(testthat)
library(ri)

####################
# Test 1. Create test matrix 1 per the PDF document.
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

# Here we set the internal arguments to createProbMatrix to ease in debugging.
raw_assignments = testmat1
byrow = F
prob_matrix = createProbMatrix(raw_assignments = raw_assignments, byrow = byrow)
prob_matrix

set.seed(4976401)
# Create sample outcome vector.
outcome = rnorm(testmat1_n)
outcome

# Choose an assignment vector for one of the permutations.
rand_column = sample(ncol(testmat1), 1)
assignment = testmat1[, rand_column]
assignment
# Make a copy equal to the internal argument of htestimate, to help with debugging.
raw_assignment = assignment

# Compare assignment 1 to assignment 2.
contrasts = c(1, -1, 0)
approx = "youngs"
result = htestimate(outcome, assignment, contrasts, prob_matrix)
# This result is an estimate of 0.93 and p-value of 0.545. Actually, now we're getting 0.607 as the p-value??
result

# This is what the variance weights are for the standard error estimate.
contrasts %*% t(contrasts)

# Test 1b - what is we change the assignments to be 0, 1, 2 rather than 1, 2, 3?
testmat1
testmat1b = testmat1 - 1
testmat1b
prob_matrix = createProbMatrix(testmat1)
prob_matrix2 = createProbMatrix(testmat1b)
prob_matrix2

# This should return true, meaning that each cell is equal. Dimnames can be different.
all.equal(prob_matrix, prob_matrix2)

# How does this compare to the ri package results?
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

####################
# Test 2. RI package example, but without blocking or clustering.
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
table(Z)
#cluster <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
#block <- c(rep(1,4),rep(2,6),rep(3,8))
#block

# Generates 10,000 samples by default.
perms = genperms(Z, maxiter=10000)
dim(perms)

# Look at the first 10 random permutations
perms[,1:10]
# Look at the assignment probability for those first 10 permutations.
colMeans(perms)[1:10]

probs <- genprob(perms) # probability of treatment
probs
ri_ate = estate(y,Z,prob=probs) # estimate the ATE
ri_ate

# Convert from 0/1 to 1/2 assignment indictators, for compatability.
z_ck = Z + 1
table(z_ck)
perms_ck = perms + 1
dim(perms_ck)
# Double-check the first 10 permutations.
perms_ck[, 1:10]
prob_matrix = createProbMatrix(perms_ck)
# This second version keeps the original 0/1 coding.
prob_matrix2 = createProbMatrix(perms)
# We should get the same probabilities between the two.
all.equal(prob_matrix, prob_matrix2)
# Examine the 1x1 matrix in the upper left.
prob_matrix[1:18, 1:18]

# Examine the htestimate results.
outcome = y
raw_assignment = z_ck
contrasts = c(-1, 1)
table(raw_assignment)
# This version manually converts the assignment levels to natural numbers, i.e. 1 and 2.
result = htestimate(y, z_ck, contrasts = c(-1, 1), prob_matrix = prob_matrix)
result
# This version keeps the original 0, 1 assignment levels.
result2 = htestimate(y, Z, contrasts = c(-1, 1), prob_matrix = prob_matrix2)
result2

# These should be true.
unlist(result) == unlist(result2)
all.equal(unlist(result), unlist(result2))

# ERROR: we are getting slightly different estimate results. What's the deal??
result$estimate == ri_ate

# Try with 100k perms
perms_100k = genperms(Z, maxiter=100000)
# It stops at 43,758 because that is all of the permutations.
choose(18, 10)
prob_matrix_100k = createProbMatrix(perms_100k)
result = htestimate(y, z_ck, contrasts = c(-1, 1), prob_matrix = prob_matrix_100k)
result

# Retry with RI.
probs_100k = genprob(perms_100k) # probability of treatment
probs_100k
ri_ate = estate(y,Z,prob=probs_100k) # estimate the ATE
ri_ate

# Confirm equality within epsilon.
abs(result$estimate - ri_ate) <= 0.0000000001
# Even more precise:
abs(result$estimate - ri_ate) <= .Machine$double.eps*2


####################
# Test 3. Compare to the RI package example, this time with clustering and blocking:
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
cluster <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
block <- c(rep(1,4),rep(2,6),rep(3,8))
block
probs <- genprobexact(Z,blockvar=block, clustvar=cluster) # probability of treatment
probs
ate <- estate(y,Z,prob=probs) # estimate the ATE
ate

#################
# TODO: Check that the Totals estimation is correct.

