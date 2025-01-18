#------------------------------------------------------------------------
#  Tests for staggered ife
#------------------------------------------------------------------------
library(tidyr)
library(ptetools)
library(BMisc)
library(dplyr)

#------------------------------------------------------------------------
#  Basic functionality
#------------------------------------------------------------------------
test_that("basic functionality", {
  sp <- reset.sim()
  sp$nife_actual <- -1
  data <- gen_ife_data(sp)

  res_lev <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = -1
  )

  sp$nife_actual <- 0
  data <- gen_ife_data(sp)
  res0 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 0
  )

  res1 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 1
  )

  sp$nife_actual <- 1
  data <- gen_ife_data(sp)
  res2 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 2
  )

  expect_equal(res_lev$overall_att$overall.att, 0, tolerance = 0.5)
  dyn_idx <- res0$event_study$egt == 0
  expect_equal(res0$event_study$att.egt[dyn_idx], 0, tolerance = 0.5)
  expect_equal(res1$overall_att$overall.att, 0, tolerance = 0.5)
  dyn_idx <- res2$event_study$egt == 0
  expect_equal(res2$event_study$att.egt[dyn_idx], 0, tolerance = 0.5)
})

test_that("empirical bootstrap", {
  sp <- reset.sim()

  sp$nife_actual <- -1
  data <- gen_ife_data(sp)
  res_lev <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = -1,
    boot_type = "empirical"
  )

  sp$nife_actual <- 0
  data <- gen_ife_data(sp)
  res0 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 0,
    boot_type = "empirical"
  )

  sp$nife_actual <- 1
  data <- gen_ife_data(sp)
  res1 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 1,
    boot_type = "empirical"
  )

  sp$nife_actual <- 2
  data <- gen_ife_data(sp)
  res2 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 2,
    boot_type = "empirical"
  )

  expect_equal(res_lev$overall_att$overall.att, 0, tolerance = 0.5)
  dyn_idx <- res0$event_study$egt == 0
  expect_equal(res0$event_study$att.egt[dyn_idx], 0, tolerance = 0.5)
  expect_equal(res1$overall_att$overall.att, 0, tolerance = 0.5)
  dyn_idx <- res2$event_study$egt == 0
  expect_equal(res2$event_study$att.egt[dyn_idx], 0, tolerance = 0.5)
})

test_that("unbalanced groups", {
  sp <- reset.sim()
  data <- gen_ife_data(sp)

  data <- subset(data, G != 6)

  res1 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 1
  )

  expect_equal(res1$overall_att$overall.att, 0, tolerance = 0.5)
})

#------------------------------------------------------------------------
#  covariates
#------------------------------------------------------------------------
test_that("include covariates", {
  skip_if(TRUE, "covariates not yet supported")
  sp <- reset.sim()
  data <- gen_ife_data(sp)

  res1 <- staggered_ife2(
    yname = "Y",
    gname = "G",
    tname = "period",
    idname = "id",
    data = data,
    nife = 1,
    xformla = ~X
  )
})

#------------------------------------------------------------------------
#  repeated cross sections
#------------------------------------------------------------------------
test_that("repeated cross sections", {
  skip_if(TRUE, "repeated cross sections not currently supported")
})

#------------------------------------------------------------------------
#  anticipation
#------------------------------------------------------------------------
test_that("anticipation", {
  skip_if(TRUE, "anticipation not currently supported")
})

#------------------------------------------------------------------------
#  unbalanced panel
#------------------------------------------------------------------------
test_that("unbalanced panel", {
  skip_if(TRUE, "unbalanced panel is not currently supported")
})
