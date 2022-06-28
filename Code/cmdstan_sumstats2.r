# Script to generate sumstats and mcmc matrix similar to rstan,
#  after fitting a model using cmdstan 
#  NOTE: following lines show sample of cmdstan fitting code
# ------------------------------------------------------
# require(parallel)
# require(cmdstanr)
# require(posterior)
# nburnin = 500                    # number warm-up (burn-in) samples 
# nsamples = 10000                 # desired total number samples
# fitmodel = c("filename.stan")    # name of file with stan code
# parms = c("Par1","Par2")         # vector of parameters to save
# stan.data <- list(N=N,X=X,Y=Y)   # list of input data variables
# cores = detectCores()
# ncore = min(40,cores-1)
# Niter = round(nsamples/ncore)
# mod <- cmdstan_model(fitmodel)   # compiles model (if necessary) 
# suppressMessages(                # Suppress messages/warnings (if desired)
#   suppressWarnings ( 
#     fit <- mod$sample(
#       data = stan.data,
#       seed = 1234,
#       chains = ncore,
#       parallel_chains = ncore,
#       refresh = 100,
#       iter_warmup = nburnin,
#       iter_sampling = Niter,
#       max_treedepth = 12,
#       adapt_delta = 0.8
#     )
#   )
# )
# ------------------------------------------------------
sumstats2 = as.data.frame(fit2$summary(variables = parms2))
row.names(sumstats2) = sumstats2$variable; sumstats2 = sumstats2[,-1] 
tmp = as.data.frame(fit2$summary(variables = parms2, mcse = mcse_mean, 
                                ~quantile(.x, probs = c(0.025, 0.975))))
sumstats2$mcse = tmp$mcse 
sumstats2$q2.5 = tmp$`2.5%` 
sumstats2$q97.5 = tmp$`97.5%`
sumstats2$q50 = sumstats2$median 
sumstats2$N_eff = sumstats2$ess_bulk
col_order = c("mean","mcse","sd","q2.5","q5","q50","q95","q97.5","N_eff","rhat")
sumstats2 = sumstats2[, col_order]
mcmc2 = as_draws_matrix(fit2$draws(variables = parms2))
vn2 = colnames(mcmc2); vns2 = row.names(sumstats2)
Nsims2 = nrow(mcmc2)
rm(tmp,col_order)   # remove temporary variables
