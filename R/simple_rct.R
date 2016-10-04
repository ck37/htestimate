
# Create a randomization distribution for a simple RCT with a certain treatment probability.
simple_rct = function(sample_size, treatment_prob, assignments = c(0, 1)) {
  counts = c(control=NA, treat=round(sample_size * treatment_prob))
  counts[1] = n - counts[2]
  counts
  # Create the box.
  box = rep(assignments, times = counts)
  table(box)

  num_perms = choose(sample_size, counts[1])
  cat("Total permutations:", prettyNum(num_perms, big.mark=","), "\n")
  # Generate the permutations.
  # TODO: suppress message about too many permutations.
  capture.output({
    perms = genperms(box)
  })
  if (ncol(perms) < num_perms) {
    cat("Switched to approximation. Sampled permutations:", prettyNum(ncol(perms), big.mark=",")
        , "\n")
  }
  return(perms)
}
