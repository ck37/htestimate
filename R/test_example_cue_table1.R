# Example from Class of Unbiased Estimators, Table 1, p. 147 - blocking and clustering.
# This function generates the assignment permutations for use in htestimate.
# This is a function so that we can re-use the cost in testing htestimate() and createProbMatrix().
test_example_cue_table1 = function() {
  y = c(1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0)
  cluster = c(1,1,2,2,3,4,5,5,5,6,6,7,7,8,9,10)
  block = c(rep(1, 6), rep(2, 10))
  x = c(4,0,4,1,4,2,4,1,2,5,4,1,4,2,2,3)

  # Create an example random asssignment vector for use in ri::genperms
  z = c(1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0)

  # Use ri package to easily generate all of the permutations.
  perms = ri::genperms(z, blockvar=block, clustvar=cluster)

  #data = data.frame(cbind(block, cluster, x, y))
  data = data.frame(block, cluster, x, y, z)

  results = list(assign_perms=perms, data=data)
  return(results)
}
