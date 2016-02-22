library(crank) # We use the permute() function.
library(testthat)
library(ri)

# Re-load htestimate.R in case we haven't already run it.
source("R/htestimate.R")

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
context("htestimate - Simple test case")

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
context("htestimate - RI package example simplified")

y = c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
table(Z)

# We skip the cluster and block data now because that will be used in a later test.

# Generates 10,000 samples by default.
# Capture warning message output from genperms.
zz = capture.output({ perms = genperms(Z, maxiter=10000) })
dim(perms)

# Look at the first 10 random permutations
perms[, 1:10]
# Look at the assignment probability for those first 10 permutations.
colMeans(perms)[1:10]

probs = genprob(perms) # probability of treatment
probs
ri_ate = estate(y, Z, prob=probs) # estimate the ATE
ri_ate

# Convert from 0/1 to 1/2 assignment indictators, for compatability.
z_ck = Z + 1
table(z_ck)
perms_ck = perms + 1
dim(perms_ck)
# Double-check the first 10 permutations.
perms_ck[, 1:10]

# This version does not need assignment to be rescaled.
prob_matrix = createProbMatrix(perms_ck)

# This second version keeps the original 0/1 coding and needs the assignment levels to be converted to 1/2.
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

# Check for symmetry in the covariance matrix.
isSymmetric(result$covariances)

# We are getting slightly different estimate results, but this is due to the # of permutations being small.
result$estimate == ri_ate

# This version keeps the original 0, 1 assignment levels.
result2 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2)
result2

# Check for symmetry in the covariance matrix.
isSymmetric(result2$covariances)

# These should be true.
unlist(result) == unlist(result2)
all.equal(unlist(result), unlist(result2))

# Constant treatment effects test.
approx = "constant effects"
result3 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2, approx = approx)
result3

# Check for symmetry in the covariance matrix.
isSymmetric(result3$covariances)

# Sharp null test.
approx = "sharp null"
result4 = htestimate(y, Z, contrasts = contrasts, prob_matrix = prob_matrix2, approx = approx)
result4

# Check for symmetry in the covariance matrix.
isSymmetric(result4$covariances)



# Try with 100k perms
perms_100k = genperms(Z, maxiter=100000)
# It stops at 43,758 because that is all of the permutations.
choose(18, 10)
prob_matrix_100k = createProbMatrix(perms_100k)
result = htestimate(y, z_ck, contrasts = c(-1, 1), prob_matrix = prob_matrix_100k)
result
# Check for symmetry in the covariance matrix.
isSymmetric(result$covariances)

# Retry with RI.
probs_100k = genprob(perms_100k) # probability of treatment
probs_100k
ri_ate = estate(y, Z, prob=probs_100k) # estimate the ATE
ri_ate

# Confirm equality within twice epsilon.
test_that("htestimate - Ex2: replicate results from RI package (no clustering or blocking).", {
  expect_lte(abs(result$estimate - ri_ate), .Machine$double.eps * 2)
})


####################
# Test 3. Compare to the RI package, with example data from Joel.
Z = c(0,0,1,1,0,0,1,1,0,0)
y = c(2,2,1,0,2,2,2,0,0,0)
block = c(rep(1,4), rep(2,6))
N = 10  # Dealing with cluster totals so N greater than M
perms = genperms(Z, blockvar=block)
dim(perms)

probs = genprobexact(Z, blockvar=block)
probs

prob_matrix = createProbMatrix(perms)
dim(prob_matrix)

# Test createProbMatrix using matrix manipulations.

# Create a stacked matrix of indicators for each assignment level.
# First indicators for control, then indicators for treatment.
tmp = rbind(1 - perms, perms)
tmpprob = tmp %*% t(tmp) / ncol(perms)
# Should get 400 TRUEs.
table(tmpprob == prob_matrix)

# Get sharp null sd using matrix operations.
ylong = c(-y ,y)
# diag(diag(x)) extracts the diagonal entries and then converts back to square matrix where
# non-diagonal entries are 0.
yexp = solve(diag(diag(tmpprob))) %*% ylong
sdest_sn = sqrt(sum(t(yexp) %*% tmpprob %*% yexp) / N^2)
sdest_sn

# Compare to ri package results under the sharp null.
Ys = genouts(y, Z, ate=0) # generate potential outcomes under sharp null of no effect
ate = estate(y, Z, prob=probs)
distout = gendist(Ys, perms, prob=probs, HT=T) # generate sampling dist. under sharp null
ri_result = dispdist(distout, ate)

# Check if they are within 2 * epsilon, which should return true.
abs(sdest_sn - ri_result$sd) <= .Machine$double.eps * 2
# This shows that the matrix manipulation and ri software generate the correct result.

# Now compare to htestimate.
ht = htestimate(y, Z, contrasts = c(-1, 1), prob_matrix, approx = "sharp null")
ht

####################
# Test 4. Compare to the RI package example, this time with clustering and blocking:
y = c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
cluster = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
block = c(rep(1, 4), rep(2, 6), rep(3, 8))
block
probs = genprobexact(Z, blockvar=block, clustvar=cluster) # probability of treatment
probs
ate = estate(y, Z, prob=probs) # estimate the ATE
ate

perms = genperms(Z, blockvar=block, clustvar=cluster)
dim(perms)

# Generate potential outcomes under tau = ate_hat
Ys = genouts(y, Z, ate=ate)
distout = gendist(Ys, perms, prob=probs)
dispdist(distout, ate)

#####
# sharp null hypothesis.
Ys = genouts(y, Z, ate=0)
distout = gendist(Ys, perms, prob=probs)
dispdist(distout, ate)

### Compare to htestimate results.
prob_matrix = createProbMatrix(perms)
htestimate(y, Z, contrasts = c(-1, 1), prob_matrix)
htestimate(y, Z, contrasts = c(-1, 1), prob_matrix, approx = "constant effects")
htestimate(y, Z, contrasts = c(-1, 1), prob_matrix, approx = "sharp null")

#################
# TODO: Check that the Totals estimation is correct.

#################
# Test 5. CUE example (table 1, p. 147), WITH blocking AND clustering.
context("htestimate - CUE table 1")
ex = test_example_cue_table1()

# This should give us a 32x32 matrix of assignment probabilities.
prob_matrix = createProbMatrix(ex$assign_perms)
dim(prob_matrix)

# We set these two variables for debugging within the function.
contrasts = c(-1, 1)
outcome = ex$y

# Loop over each possible assignment permutation and calculate a separate result.
results = list()
for (perm_i in 1:ncol(ex$assign_perms)) {
  assignment = ex$assign_perms[, perm_i]
  # Calculate the HT estimate of the treatment effect.
  result = htestimate(ex$data$y, assignment, contrasts=contrasts, prob_matrix)
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

#head(results)

# ERROR: we are getting NaNs for some of the std errors, because the covariances are much larger than the variances.
# So we may have an error in either the covariance or the variance term calculations, or perhaps the weight calculation?

# See e.g. assignment = assign_perms[, 78]
assignment = ex$assign_perms[, 78]
result = htestimate(ex$data$y, assignment, contrasts=contrasts, prob_matrix)
result
# Check for symmetry in the covariance matrix.
isSymmetric(result$covariances)

# What if we use a different covariance approximation?
result = htestimate(ex$data$y, assignment, contrasts=contrasts, prob_matrix, approx="constant effects")
result
# Check for symmetry in the covariance matrix.
isSymmetric(result$covariances)

# Compare to sharp null results.
htestimate(ex$data$y, assignment, contrasts=contrasts, prob_matrix, approx="sharp null")


## Replicating table 2.

# Expected value of estimate, should be 0 per table 2, p. 149 of CUE.
estimates = sapply(results, FUN=function(x){ x$estimate })
mean(estimates)
# Standard error of delta from the paper (first SE row in table 2, under "HT" column).
sqrt(sum((estimates - mean(estimates))^2)/length(estimates))
# This does match the table's result.

# This is the actual variance of the estimated ATE.
sum((estimates - mean(estimates))^2)/length(estimates)
errors = sapply(results, FUN=function(x){ x$std_err })

# This should be 0, but we are getting 9 NaNs due to negative terms in the covariance matrix.
sum(is.na(errors))

# This version gives us the expectation of the estimated variance.
mean(errors^2)
sum((errors - mean(errors))^2)/length(errors)

# TODO: clean up test.

context("htestimate - rerandomization 1")

# Simulate a sample dataset with 1 covariate and 1 outcome.
n = 20
set.seed(1)
# x1 is our observed covariate.
x1 = runif(n)
summary(x1)
# U1 is an unobserved variable.
u1 = rexp(n)
# E is our random error.
e = rnorm(n)
# t is our treatment effect.
t = 2
y0 = x1 + 2 * u1 + e
# Constant treatment effect:
y1 = y0 + t

data = data.frame(x1, u1, y0, y1)

# How many units to allocation to each assignment level.

# Generate randomization distribution
perms = simple_rct(20, 0.5)

# Analyze balance based on the x1 variable.
rand_analyzed = analyze_randomizations(perms, data.frame(data[, "x1"]))

rownames(rand_analyzed)

rand_f0.2 = restrict_randomizations(rand_analyzed, 0.2)
rand_f0.9 = restrict_randomizations(rand_analyzed, 0.9)


dim(rand_analyzed)
rownames(rand_analyzed)
summary(rand_analyzed["f_p", ])
summary(rand_analyzed["r_sqr", ])

table(rand_analyzed["keep", ])
prop.table(table(rand_analyzed["keep", ]))
# Look at the distribution of r_squared for keep vs discard.
tapply(rand_analyzed["f_p", ], rand_analyzed["keep", ], FUN=summary)

# Select which randomization permutations we want to keep.
perms_restricted = perms[, as.logical(rand_analyzed["keep", ])]
dim(perms_restricted)

# Choose one permutation to be "observed".
prob_matrix = createProbMatrix(perms_restricted)
prob_matrix

prob_matrix = createProbMatrix(perms)
prob_matrix
# The diagonal is in theory 0.5, but we see some finite sample variance.
mean(diag(prob_matrix))
sd(diag(prob_matrix))

# Ratio of the variances of two designs.
#
context("joel - test case 1")

Z<-c(0,0,1,1,0,0,1,1,0,0)
y<-c(2,2,1,0,2,2,2,0,0,0)
block<-c(rep(1,4), rep(2,6))
N<-10 #dealing with cluster totals so N greater than M
perms<-genperms(Z, blockvar=block)
probs<-genprobexact(Z,blockvar=block)

#test Chris' software
prob_matrix<-createProbMatrix(perms)

#test createProbMatrix
tmp<-rbind(1-perms,perms)
tmpprob<-tmp%*%t(tmp)/ncol(perms)
table(tmpprob==prob_matrix)

#get sharp null sd using matrix operations
ylong<-c(-y,y)
# y_expanded = y over the probabilities (P^-1 * Y).
yexp = solve(diag(diag(tmpprob))) %*% ylong
sdest_sn<-sqrt(sum(t(yexp) %*% tmpprob %*% yexp)/N^2)
sdest_sn


#exactly same as Peter's software under sharp null
Ys <- genouts(y,Z,ate=0) # generate potential outcomes under sharp null of no effect
ate<-estate(y,Z, prob=probs)
ate
distout <- gendist(Ys,perms, prob=probs, HT=T) # generate sampling dist. under sharp null
dispdist(distout, ate)
# Compare to htestimate:
htestimate(y, Z, contrasts = c(-1, 1), prob_matrix=prob_matrix, approx="sharp null",  totals=F)
