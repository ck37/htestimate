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
assignments = testmat1
byrow = F
prob_matrix = createProbMatrix(assignments = assignments, byrow = byrow)
prob_matrix

set.seed(4976401)
# Create sample outcome vector.
outcome = rnorm(testmat1_n)
outcome

# Choose an assignment vector for one of the permutations.
rand_column = sample(ncol(testmat1), 1)
assignment = testmat1[, rand_column]
assignment

# Compare assignment 1 to assignment 2.
contrasts = c(1, -1, 0)
approx = "youngs"
result = htestimate(outcome, assignment, contrasts, prob_matrix, approx=approx)
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
context("RI package example simplified")

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
assignment = z_ck
contrasts = c(-1, 1)
totals = F
table(assignment)
# This version manually converts the assignment levels to natural numbers, i.e. 1 and 2.
result = htestimate(y, z_ck, contrasts = contrasts, prob_matrix = prob_matrix)
result
# This version keeps the original 0, 1 assignment levels.
result2 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2)
result2

# These should be true.
unlist(result) == unlist(result2)
all.equal(unlist(result), unlist(result2))

# Constant treatment effects test.
approx = "constant effects"
result3 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2, approx = approx)
result3

# Sharp null test.
approx = "sharp null"
result4 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2, approx = approx)
result4

# We are getting slightly different estimate results, but this is due to the # of permutations being small.
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

# TODO: add in htestimate version.

#################
# TODO: Check that the Totals estimation is correct.

#################
# Test 4. CUA example (table 1, p. 147), WITH blocking AND clustering.
y = c(1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0)
cluster = c(1,1,2,2,3,4,5,5,5,6,6,7,7,8,9,10)
block = c(rep(1, 6), rep(2, 10))
x = c(4,0,4,1,4,2,4,1,2,5,4,1,4,2,2,3)

data = data.frame(cbind(block, cluster, x, y))
data = data.frame(block, cluster, x, y)

# Collapse to the cluster level before performing randomization.
library(dplyr)
data_clusters = distinct(data, cluster) %>% select(block, cluster)
data_clusters

# Generate many possible blocked random assignments; each is a column.
library(randomizr)
set.seed(1)
# Assign 2 units to treatment per block (50% chance in 1st block, 33% chance in 2nd block)
block_m = rbind(c(2, 2),
                 c(4, 2))
reps = replicate(900, block_ra(data_clusters$block, block_m = block_m, condition_names = c("control", "treatment")))
# Set margin=2 so that matrix is unique by column (assignment permutation) rather than row.
cluster_perms = unique(reps, MARGIN=2)
# We expect and do get 90 unique possible assignments.
dim(cluster_perms)

# Review the Pi_1i's by cluster - probability of being assigned treatment.
# These should be 0.5 for the first 4 clusters then 0.333 for the last 6.
apply(cluster_perms, MARGIN=1, FUN=function(x){ mean(x == "treatment") })

# Expand the cluster assignment to the unit level.
assign_perms = apply(cluster_perms, MARGIN=2, FUN=function(assignment) {
  # Create a dataframe with the cluster ID and the cluster-level assignment.
  assign_df = data.frame(cluster=unique(data$cluster), assignment)
  # Merge cluster assignment back to the unit-level data.
  unit_assignment = merge(data, assign_df, by="cluster")
  # We convert to numeric rather than a factor to see if this can be fix the NAs in covariance bug.
  # But it seems not to be helping right now.
  as.numeric(as.factor(unit_assignment$assignment))
})
# We should have a 16x90 matrix of assignment permutations at the unit level.
dim(assign_perms)

# Confirm that our treatment assignment probabilities remain 0.5 and then 1/3.
apply(assign_perms, MARGIN=1, FUN=function(x){ mean(x == "treatment") })
# We need to use this version once we have converted to numeric rather than string values.
apply(assign_perms, MARGIN=1, FUN=function(x){ mean(x == 2) })

# This should give us a 32x32 matrix of assignment probabilities.
prob_matrix = createProbMatrix(assign_perms)
dim(prob_matrix)

# We set these two variables for debugging within the function.
contrasts = c(-1, 1)
outcome = y

# Loop over each possible assignment permutation.
results = list()
for (perm_i in 1:ncol(assign_perms)) {
  assignment = assign_perms[, perm_i]
  # Calculate the HT estimate of the treatment effect.
  result = htestimate(y, assignment, contrasts=contrasts, prob_matrix)
  if (is.nan(result$std_err)) {
    cat("Error in permutation", perm_i, "\n")
    cat("Assignment:")
    print(assignment)
    cat("Result:\n")
    print(result)
    cat("\n")
  }
  results[[perm_i]] = result
}

head(results)

# ERROR: we are getting NaNs for some of the std errors, presumably because the variance is negative :/.


## Replicating table 2.

# Expected value of estimate, should be 0 per table 2, p. 149 of CUE.
estimates = sapply(results, FUN=function(x){ x$estimate })
mean(estimates)
# Standard error of delta from the paper (first SE row in table 2).
sqrt(sum((estimates - mean(estimates))^2)/length(estimates))
# This does match the table's result.

# This is the actual variance of the estimated ATE.
sum((estimates - mean(estimates))^2)/length(estimates)
errors = sapply(results, FUN=function(x){x$std_err})
# This version gives us the expectation of the estimated variance.
mean(errors^2, na.rm=T)
sum((errors - mean(errors))^2)/length(errors)
