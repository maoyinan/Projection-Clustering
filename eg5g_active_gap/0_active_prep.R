# manually create gaps in time series to test method

library(tidyverse)
filter <- dplyr::filter

#' ### Import dataset
load('../eg5_active/dat')

set.seed(123)
idx <- sample(nrow(dat), 650, replace=F)
idx

dat <- dat[-idx,]

save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')


# linear mixed model ------------------------------------
source('../1_main.R')

out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg7Fixed')

save(df_of_draws,file='D:/df_of_drawsEg7')
