% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ife.R, R/imports.R
\docType{package}
\name{ife}
\alias{ife}
\title{ife}
\usage{
ife(
  yname,
  gname,
  tname,
  idname,
  data,
  nife,
  xformla = ~1,
  zformla,
  ret_ife_regs = TRUE,
  anticipation = 0,
  cband = TRUE,
  alp = 0.05,
  boot_type = "multiplier",
  biters = 100,
  cl = 1
)
}
\arguments{
\item{nife}{the number of interactive fixed effects to include in the model}
}
\value{
\code{pte::pte_results} object
}
\description{
Compute treatment effects in interactive fixed effects models
with a small number of time periods

The \code{ife} package contains tools for estimating treatment
effects when untreated potential outcomes are generated by an interactive
fixed effects model.
}