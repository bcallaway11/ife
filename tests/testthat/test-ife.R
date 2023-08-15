
# clear console: ctrl+L
#testthat::test_dir("tests/testthat")
#testthat::test_file("tests/testthat/test-ife.R")

# Download some packages
# devtools::install_github("bcallaway11/pte")
# install.packages("ivreg")


# Load functions
fldr <- "C:/Users/yjeon/DOCUME~1/MOBAXT~1/home/ife/R/"
sapply(paste0(fldr,list.files(fldr)), source)

#------------------------------------------------------------------------
#  Tests for ife
#-------------------------------------------------------------------------------
library(tidyr)
library(BMisc)
library(dplyr)
library(pte)

#-------------------------------------------------------------------------------
#  Basic functionality
#-------------------------------------------------------------------------------
# basic functionality ==========================================================
test_that("basic functionality", {
  sp <- reset.sim2()

## number of actual ife = 0 ----------------------------------------------------
  sp$nife_actual <- 0
  data <- gen_ife_data2(sp)

  # number of ife (tried) = 0
  res00 <- ife(yname="Y",
             gname="G",
             tname="period",
             idname="id",
             data=data,
             nife=0,
             xformla=~1,
             zformla=~1)
  
  expect_equal(res00$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res00$event_study$egt == 0
  expect_equal(res00$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  # number of ife (tried) = 1 & Z1
  res01_1 <- ife(yname="Y",
              gname="G",
              tname="period",
              idname="id",
              data=data,
              nife=1,
              xformla=~1,
              zformla=~1+Z1)
  
  expect_equal(res01_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res01_1$event_study$egt == 0
  expect_equal(res01_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  # number of ife (tried) = 1 & Z2
  res01_2 <- ife(yname="Y",
               gname="G",
               tname="period",
               idname="id",
               data=data,
               nife=1,
               xformla=~1,
               zformla=~1+Z2)
  
  expect_equal(res01_2$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res01_2$event_study$egt == 0
  expect_equal(res01_2$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  

  # number of ife (tried) = 1 & Z3
  res01_3 <- ife(yname="Y",
               gname="G",
               tname="period",
               idname="id",
               data=data,
               nife=1,
               xformla=~1,
               zformla=~1+Z3)
  
  expect_equal(res01_3$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res01_3$event_study$egt == 0
  expect_equal(res01_3$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  

## number of actual ife = 1 ----------------------------------------------------
  sp$nife_actual <- 1
  data <- gen_ife_data2(sp)
  
  # number of ife (tried) = 1 & Z1
  res11_1 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=1,
                 xformla=~1,
                 zformla=~1+Z1)
  
  expect_equal(res11_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res11_1$event_study$egt == 0
  expect_equal(res11_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  
## number of actual ife = 2 ----------------------------------------------------
  sp$nife_actual <- 2
  data <- gen_ife_data2(sp)
  
  # number of ife (tried) = 2 & Z1, Z2
  res22_12 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=2,
                 xformla=~1,
                 zformla=~1+Z1+Z2)
  
  expect_equal(res22_12$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res22_12$event_study$egt == 0
  expect_equal(res22_12$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  
## number of actual ife = 3 <- this fails most of the time...???!!! ----------------------------------------------------
  sp$nife_actual <- 3
  data <- gen_ife_data2(sp)
  
  # number of ife (tried) = 3 & Z1, Z2, Z3
  res33_123 <- ife(yname="Y",
                  gname="G",
                  tname="period",
                  idname="id",
                  data=data,
                  nife=3,
                  xformla=~1,
                  zformla=~1+Z1+Z2+Z3)
  
  expect_equal(res33_123$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res33_123$event_study$egt == 0
  expect_equal(res33_123$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  
}) # basic functionality



# unbalanced groups ============================================================
test_that("unbalanced groups", {
  
  sp <- reset.sim2()
  
## number of actual ife = 0 ----------------------------------------------------
  sp$nife_actual <- 0
  data <- gen_ife_data2(sp)
  data <- subset(data, G != 6)
  
  # number of ife (tried) = 0
  res00 <- ife(yname="Y",
               gname="G",
               tname="period",
               idname="id",
               data=data,
               nife=0,
               xformla=~1,
               zformla=~1)
  
  expect_equal(res00$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res00$event_study$egt == 0
  expect_equal(res00$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  # number of ife (tried) = 1 & Z1
  res01_1 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=1,
                 xformla=~1,
                 zformla=~1+Z1)
  
  expect_equal(res01_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res01_1$event_study$egt == 0
  expect_equal(res01_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  

  
## number of actual ife = 1 ----------------------------------------------------
  sp$nife_actual <- 1
  data <- gen_ife_data2(sp)
  data <- subset(data, G != 6)
  
  # number of ife (tried) = 1 & Z1
  res11_1 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=1,
                 xformla=~1,
                 zformla=~1+Z1)
  
  expect_equal(res11_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res11_1$event_study$egt == 0
  expect_equal(res11_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  
}) # unbalanced groups



# empirical bootstrap
test_that("empirical bootstrap", {
  sp <- reset.sim2()

## number of actual ife = 0 ----------------------------------------------------
  sp$nife_actual <- 0
  data <- gen_ife_data2(sp)

  # number of ife (tried) = 0
  res00 <- ife(yname="Y",
               gname="G",
               tname="period",
               idname="id",
               data=data,
               nife=0,
               xformla=~1,
               zformla=~1,
               boot_type="empirical")

  expect_equal(res00$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res00$event_study$egt == 0
  expect_equal(res00$event_study$att.egt[dyn_idx], 0, tolerance=0.5)

}) # empirical bootstrap
  
#------------------------------------------------------------------------
#  covariates (X)
#------------------------------------------------------------------------
# covariates (X)
test_that("include covariates", {
  sp <- reset.sim2()
  
## number of actual ife = 0 ----------------------------------------------------
  sp$nife_actual <- 0
  data <- gen_ife_data2(sp)
  
  # number of ife (tried) = 0
  res00 <- ife(yname="Y",
               gname="G",
               tname="period",
               idname="id",
               data=data,
               nife=0,
               xformla=~1+X1+X2,
               zformla=~1+X1+X2)
  
  expect_equal(res00$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res00$event_study$egt == 0
  expect_equal(res00$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  # number of ife (tried) = 1 & Z1
  res01_1 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=1,
                 xformla=~1+X1+X2,
                 zformla=~1+X1+X2+Z1)
  
  expect_equal(res01_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res01_1$event_study$egt == 0
  expect_equal(res01_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
  
## number of actual ife = 1 ----------------------------------------------------
  sp$nife_actual <- 1
  data <- gen_ife_data2(sp)
  
  # number of ife (tried) = 1 & Z1
  res11_1 <- ife(yname="Y",
                 gname="G",
                 tname="period",
                 idname="id",
                 data=data,
                 nife=1,
                 xformla=~1+X1+X2,
                 zformla=~1+X1+X2+Z1)
  
  expect_equal(res11_1$overall_att$overall.att, 0, tolerance=0.5) # Differences smaller than tolerance are not reported.
  dyn_idx <- res11_1$event_study$egt == 0
  expect_equal(res11_1$event_study$att.egt[dyn_idx], 0, tolerance=0.5)
  
}) # covariates (X)

