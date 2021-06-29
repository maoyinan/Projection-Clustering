library(rstan)
library(rstanarm)
library(parallel)
library(rstudioapi)
library(grid)
library(tidyverse)
library(kableExtra)
library(foreach)
library(doParallel)
options(mc.cores = parallel::detectCores())

source('../2_postFunction.R')
source('../3_supportFunctions.R')
source('../4_plotCluster.R')
source('../5_plotKmeans_fitted.R')
source('../5_plotKmeans.R')
source('../6_plotFit.R')

# Import dataset
load('dat')

# Fit Bayesian mixed linear regression
# Use lmer to get prior distribution.

y_var <- 'Record'
x_vars <- paste0('Z',1:nBasis)
id_var <- 'ID'
# fml <- sprintf('%s~  0 + (0 + %s|%s)', y_var, paste(x_vars, collapse=' + '),id_var)
# fml <- sprintf('%s~  (0 + %s|%s)', y_var, paste(x_vars, collapse=' + '),id_var)
# fml <- as.formula(fml)
# fml
# temp1 <- stan_lmer(
#   fml,
#   prior_covariance = decov(1),
#   prior_aux= exponential(rate=1),
#   dat, chains=4, seed=1)
# prior_summary(temp1)
# dat$y_hat <- fitted(temp1)

# ------
#
#   Auxiliary (sigma)
# Specified prior:
#   ~ exponential(rate = 1)
# Adjusted prior:
#   ~ exponential(rate = 1)
#
# Covariance
# ~ decov(reg. = 1, conc. = 1, shape = 1, scale = 1)
# ------
#   See help('prior_summary.stanreg') for more details


# MCMC algorithm to approximate Bayesian mixed linear regression
dat_mcmc <- list(
  N= nrow(dat),
  y= dat$Record,
  nBasis= nBasis,
  Z= dat[,paste0('Z',seq(nBasis))],
  nSub= length(unique(dat$ID)),
  sub= as.numeric(dat$ID)
)
rt <- stanc(file="../mcmc_try2.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(fit02 <- sampling(sm, data=dat_mcmc, seed=1,cores=detectCores(), init = 0))
df_of_draws <- as.data.frame(fit02)
# save(df_of_draws,file='df_of_draws')

# Projection based clustering
# load('df_of_draws')

ls_par <- post.mean(df_of_draws, dat, nBasis)
mat_fitted <- pc.fit(ls_par, dat, ls_idxA,seed=1)$mat_fitted




