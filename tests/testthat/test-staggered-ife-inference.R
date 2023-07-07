#------------------------------------------------------------------------
#  Tests related to inference for staggered interactive fixed effects
#------------------------------------------------------------------------

library(tidyr)
library(pte)
library(BMisc)
library(dplyr)
library(pbapply)

cl <- as.numeric(readline("how many clusters would like to use for the inference tests?  "))

test_that("tests for inference", {
  mc_sims <- 100
  
  #------------------------------------------------------------------------
  #  levels
  #------------------------------------------------------------------------
  rejs <- pblapply(1:mc_sims, function(mc) {
    sp <- reset.sim()
    sp$nife_actual <- -1
    data <- gen_ife_data(sp)
    res <- staggered_ife2(yname="Y",
                          gname="G",
                          tname="period",
                          idname="id",
                          data=data,
                          nife=-1)
    tstat_overall <- res$overall_att$overall.att / res$overall_att$overall.se
    tstat_attgt <- res$att_gt$att[5] / res$att_gt$se[5]
    tstat_dyn <- res$event_study$att.egt[4] / res$event_study$se.egt[4]
    rej_overall <- 1*( abs(tstat_overall) > qnorm(.975))
    rej_attgt <- 1*( abs(tstat_attgt) > qnorm(.975))
    rej_dyn <- 1*( abs(tstat_dyn) > qnorm(.975))
    list(rej_overall=rej_overall, rej_attgt=rej_attgt, rej_dyn=rej_dyn)
  }, cl=cl)
  
  rej_overall_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_overall")))
  rej_attgt_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_attgt")))
  rej_dyn_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_dyn")))
  
  expect_equal(rej_overall_frac, 0.06, tolerance=.99) # make test fail if reject 0
  expect_equal(rej_attgt_frac, 0.06, tolerance=.99) 
  expect_equal(rej_dyn_frac, 0.06, tolerance=.99)
  
  
  #------------------------------------------------------------------------
  #  0 ife
  #------------------------------------------------------------------------
  rejs <- pblapply(1:mc_sims, function(mc) {
    sp <- reset.sim()
    sp$nife_actual <- 0
    data <- gen_ife_data(sp)
    res <- staggered_ife2(yname="Y",
                          gname="G",
                          tname="period",
                          idname="id",
                          data=data,
                          nife=0)
    tstat_overall <- res$overall_att$overall.att / res$overall_att$overall.se
    tstat_attgt <- res$att_gt$att[5] / res$att_gt$se[5]
    tstat_dyn <- res$event_study$att.egt[4] / res$event_study$se.egt[4]
    rej_overall <- 1*( abs(tstat_overall) > qnorm(.975))
    rej_attgt <- 1*( abs(tstat_attgt) > qnorm(.975))
    rej_dyn <- 1*( abs(tstat_dyn) > qnorm(.975))
    list(rej_overall=rej_overall, rej_attgt=rej_attgt, rej_dyn=rej_dyn)
  }, cl=cl)
  
  rej_overall_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_overall")))
  rej_attgt_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_attgt")))
  rej_dyn_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_dyn")))
  
  expect_equal(rej_overall_frac, 0.06, tolerance=.99) # make test fail if reject 0
  expect_equal(rej_attgt_frac, 0.06, tolerance=.99) 
  expect_equal(rej_dyn_frac, 0.06, tolerance=.99)
  
  #------------------------------------------------------------------------
  #  1 ife
  #------------------------------------------------------------------------
  rejs <- pblapply(1:mc_sims, function(mc) {
    sp <- reset.sim()
    data <- gen_ife_data(sp)
    res <- staggered_ife2(yname="Y",
                           gname="G",
                           tname="period",
                           idname="id",
                           data=data,
                           nife=1)
    tstat_overall <- res$overall_att$overall.att / res$overall_att$overall.se
    tstat_attgt <- res$att_gt$att[5] / res$att_gt$se[5]
    tstat_dyn <- res$event_study$att.egt[4] / res$event_study$se.egt[4]
    rej_overall <- 1*( abs(tstat_overall) > qnorm(.975))
    rej_attgt <- 1*( abs(tstat_attgt) > qnorm(.975))
    rej_dyn <- 1*( abs(tstat_dyn) > qnorm(.975))
    list(rej_overall=rej_overall, rej_attgt=rej_attgt, rej_dyn=rej_dyn)
  }, cl=cl)
  
  rej_overall_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_overall")))
  rej_attgt_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_attgt")))
  rej_dyn_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_dyn")))
  
  expect_equal(rej_overall_frac, 0.06, tolerance=.99) # make test fail if reject 0
  expect_equal(rej_attgt_frac, 0.06, tolerance=.99) 
  expect_equal(rej_dyn_frac, 0.06, tolerance=.99)
  
  #------------------------------------------------------------------------
  #  2 ifes
  #------------------------------------------------------------------------
  
  rejs <- pblapply(1:mc_sims, function(mc) {
    sp <- reset.sim()
    sp$nife_actual <- 2
    data <- gen_ife_data(sp)
    res <- staggered_ife2(yname="Y",
                          gname="G",
                          tname="period",
                          idname="id",
                          data=data,
                          nife=2)
    tstat_overall <- res$overall_att$overall.att / res$overall_att$overall.se
    tstat_attgt <- res$att_gt$att[5] / res$att_gt$se[5]
    tstat_dyn <- res$event_study$att.egt[4] / res$event_study$se.egt[4]
    rej_overall <- 1*( abs(tstat_overall) > qnorm(.975))
    rej_attgt <- 1*( abs(tstat_attgt) > qnorm(.975))
    rej_dyn <- 1*( abs(tstat_dyn) > qnorm(.975))
    list(rej_overall=rej_overall, rej_attgt=rej_attgt, rej_dyn=rej_dyn)
  }, cl=cl)
  
  rej_overall_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_overall")))
  rej_attgt_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_attgt")))
  rej_dyn_frac <- mean(unlist(BMisc::getListElement(rejs, "rej_dyn")))
  
  expect_equal(rej_overall_frac, 0.06, tolerance=.99) # make test fail if reject 0
  expect_equal(rej_attgt_frac, 0.06, tolerance=.99) 
  expect_equal(rej_dyn_frac, 0.06, tolerance=.99)
})