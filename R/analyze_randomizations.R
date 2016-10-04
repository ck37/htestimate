# New function to randomly permute assignments within blocks, after aggregating to cluster totals.
# Or do randomizr, blockTools, or ri (et al?) already provide this functionality?
# Examine balance in a set of randomization permutations.
analyze_randomizations = function(permutations, covariates) {
  # TODO: parallelize
  results = apply(perms, covariates, MARGIN=2, FUN = function(assignment, data) {
    # Conduct a balance test.
    #TODO: multinomial logit for > 2 treatment arms.
    reg = lm(assignment ~ ., data=data)
    reg_sum = summary(reg)
    f_stat = reg_sum$fstatistic
    # Calculate the p-value of the f-statistic, which doesn't seem otherwise extractable.
    f_p = pf(f_stat[1L], f_stat[2L], f_stat[3L], lower.tail = FALSE)
    names(f_p) = NULL
    r_sqr = reg_sum$r.squared
    #keep_randomization = as.numeric(f_p >= f_p_cutoff)
    #results = c(keep=keep_randomization, f_p=f_p, r_sqr=r_sqr)
    results = c(f_p=f_p, r_sqr=r_sqr)
    results
  })
  results
}
