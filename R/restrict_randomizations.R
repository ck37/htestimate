# Restrict our randomizations to those with balance test p-value above a threshold.
restrict_randomizations = function(analyzed_rand_perms, f_p_threshold) {
  keep = c(keep=analyzed_rand_perms["f_p", ] >= f_p_threshold)
  rand_analyzed = rbind(analyzed_rand_perms, keep)
  return(rand_analyzed)
}
